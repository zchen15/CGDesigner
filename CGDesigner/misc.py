#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# Miscellaneous functions

import json
import datetime
import time
import os
import psutil
import hashlib
import logging

import pandas as pd
import nupack

def read_config(fname):
    '''
    Reads json formatted data from a text file
    fname = input json data file
    '''
    print('reading config from', fname)
    with open(fname, 'r') as f:
        data = json.load(f)
    return data

def write_config(fname, data):
    '''
    Export data in json format to a text file
    fname = output filename
    data = dictionary of values
    '''
    print('writing config to', fname)
    with open(fname, 'w') as f:
        json.dump(data, f, indent=2)

def update_parameters(args, params=None):
    '''
    Update params dictionary with values from args dictionary
    args = dictionary of values to apply to params
    params = existing dictionary of parameters
    '''
    if params==None:
        params = {}
    if args!=None:
        # format to dictionary
        if type(args) == dict:
            val = args
        else:
            val = args.__dict__
        # update values
        for k in val.keys():
            if val[k]!=None:
                params[k] = val[k]
    return params

def get_timestamp():
    '''
    Get a time stamp in year month day hr min sec format
    '''
    return datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')

def get_memory_usage():
    '''
    Get current memory usage information
    '''
    pid = os.getpid()
    process = psutil.Process(pid)
    return process.memory_full_info().rss/1024/1024

def get_hash(x):
    '''
    Generates a hash for a given input x
    '''
    y = get_timestamp()+x
    h = hashlib.sha224(y.encode()).hexdigest()
    return h

def read_csv(fname):
    '''
    Reads csv and removes spaces from sequence
    '''
    try:
        if type(fname) == str:
            logging.info('Reading '+fname)
            df = pd.read_csv(fname)
            df['Sequence'] = [i.replace(' ','') for i in df['Sequence']]
            return df
        
        elif type(fname) == list:
            out = []
            for f in fname:
                x = read_csv(f)
                if ('filename' in x.columns)==False:
                    x['filename'] = f
                out.append(x)
            return pd.concat(out)
    except:
        logging.info('Failed to read '+str(fname))
        return []

def parse_results(infiles, outfile, data):
    '''
    Extract info from nupack output
    '''
    out = []
    for fname in infiles:
        out.append(get_dataframe(fname, data))
    out = pd.concat(out)
    out.to_csv(outfile+'.csv', index=False)

def get_dataframe(fname, field):
    '''
    Generate dataframe from nupack results
    '''
    res = nupack.DesignResult.load(fname)
    defect = res.defects.ensemble_defect
    logging.info(fname+' defect:'+str(defect))
    if field=='strands':
        # export to strands
        df = res.to_analysis.strand_table()
        df['Sequence'] = [i.replace('U','T') for i in df['Sequence']]
    elif field=='domains':
        # export domains
        df = res.to_analysis.domain_table()
        df['Sequence'] = [i.replace('U','T') for i in df['Sequence']]
    elif field=='tubes':
        # get tubes
        df = res.concentrations.table
    df['defect'] = defect
    df['filename'] = fname
    return df