#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# Class function for processing the cgRNA designs

import nupack
import numpy as np
import logging
import time
import pickle
import subprocess
import os

from .fileserver import Fileserver
from .misc import update_parameters
from .misc import get_memory_usage
from .misc import read_config

class Worker:
    params = {}
    params['s3'] = {}

    def __init__(self, args=None):
        self.params = update_parameters(args, self.params)
        # read config if available
        if 'config' in self.params:
            data = read_config(self.params['config'])
            self.params = update_parameters(data, self.params)
            del(self.params['config'])
        # set target bucket for uploading files to s3
        key = ['bucket','aws_access_key_id','aws_secret_access_key','endpoint_url']
        for k in key:
            if k in self.params:
                self.params['s3'][k] = self.params[k]
        self.server = Fileserver(self.params)

    def run(self):
        '''
        Restart designs from checkpoint files
        '''
        # find list of files
        files = self.params['infile']
        # fetch the file
        self.server.fetch_files(files)
        for i in range(len(files)):
            fname = files[i]
            if os.path.exists(fname)==False:
                try:
                    logging.info('Fetching '+fname)
                    outname = fname.split('/')[-1]
                    cmd = 'wget '+fname+' -O '+outname
                    subprocess.run(cmd.split(' '))
                    files[i] = outname
                except:
                    logging.info('Failed to fetch '+fname)
        
        files = np.array(files)
        c = np.array([('.spec' in f) for f in files])
        for f in files[~c]:
            try:
                logging.info('Running '+f)
                d = pickle.load(f)
                d.run()
            except:
                logging.info('Failed to run '+f)
            try:
                logging.info('Executing '+f)
                cmd = 'bash '+f
                subprocess.run(cmd.split(' '))
            except:
                logging.info('Failed to execute '+f)

        designs = []
        for f in files[c]:
            d = nupack.tube_design.load(f)
            name = f.split('.spec')[0]
            logging.info('adding '+name+' to design queue')
            designs.append([name, d])
        self.dispatch(designs)

    def dispatch(self, designs):
        '''
        Dispatch nupack jobs to queue
        '''
        interval = self.params['ckpt_int']
        trials = self.params['trials']
        jobs = []
        while len(designs) + len(jobs) > 0:
            jobs = self.add_jobs(jobs, designs, trials)
            jobs = self.check_jobs(jobs)
            time.sleep(interval)

    def add_jobs(self, jobs, stack, trials):
        '''
        Add jobs from stack
        jobs = currently running jobs
        stack = all jobs that need to be run
        trials = max number of jobs in jobs
        '''
        logging.info('Memory used=' + str(get_memory_usage()))
        N = trials - len(jobs)
        for i in range(N):
            if len(stack) > 0:
                n, d = stack.pop(0)
                jobs.append([n, d.launch(trials=1, interval=self.params['ckpt_int'], checkpoint=n), time.time()])
        logging.info('jobs in queue='+str(len(stack)))
        return jobs

    def check_jobs(self, jobs):
        '''
        Check on job status and save results if jobs finished
        returns list of unfinished jobs
        '''
        keep = [] 
        max_opt = self.params['max_opt']
        for t in range(len(jobs)):
            name, d, start = jobs[t]
            # print elapse time
            elapse = time.time()-start
            msg = '[name, elapse, obj]: '+str([name, elapse, d])
            logging.info(msg)
            # check if design finished or exceeded max_opt time
            if d.trials[0].future.done() or elapse > max_opt:
                logging.info('Completed '+name)
                # export results
                try:
                    d.stop()
                    res = d.final_results()[0]
                    self.checkpoint(name, res, timestamp=False, defect=True, csv=True)
                except:
                    logging.error('Failed to export results from '+name)
            else:
                keep.append(jobs[t])
        return keep

    def checkpoint(self, fname, result, timestamp=False, defect=True, csv=False):
        '''
        Save progress
        '''
        if timestamp:
            fname+= '_'+ts 
        if defect:
            fname+= '_d'+str(result.defects.ensemble_defect)
        fname+='.o'
        result.save(fname)
        self.server.upload_files(fname)
        if csv:
            fname = self.export_to_csv(fname)

    def export_to_csv(self, fname):
        '''
        Save strands to csv file
        '''
        res = nupack.DesignResult.load(fname)
        defect = res.defects.ensemble_defect
        logging.info(fname+' defect:'+str(defect))

        # export to strands 
        df = res.to_analysis.strand_table()
        df['Sequence'] = [i.replace('U','T') for i in df['Sequence']]
        
        # write to output
        fout = fname.split('.o')[0]+'_strands.csv'
        logging.info('writing strands to '+fout)
        df = df.sort_values(by=df.columns[0])
        df.to_csv(fout, index=False)
        self.server.upload_files(fout)

