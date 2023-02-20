#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# Executable for cgRNA design

from .designer import Designer
from .worker import Worker
from .fileserver import Fileserver
from .kubectrl import Kube
from .analysis import Analysis
from .oligos import Oligos
from .misc import parse_results

import numpy as np
import argparse
import os
import logging
import sys

def set_log(filename='job.log', log_level='INFO'):
    fmt = 'pid['+str(os.getpid())+'] %(asctime)s.%(msecs)03d %(levelname)s: %(message)s'
    dfmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(fmt, datefmt=dfmt)
    # print to stdout
    handle = logging.StreamHandler(sys.stdout)
    handle.setLevel(log_level)
    handle.setFormatter(formatter)
    logging.getLogger().addHandler(handle)
    # print to log file
    handle = logging.FileHandler(filename=filename)
    handle.setLevel(log_level)
    handle.setFormatter(formatter)
    logging.getLogger().addHandler(handle)
    logging.getLogger().setLevel(log_level)

def main():
    '''
    Main function
    '''
    parser = argparse.ArgumentParser(description='A python based tool for designing and analyzing conditional guide RNAs with NUPACK')
    parser.add_argument('-v', dest='verbose', action='store_true', help='print verbose info')
    parser.add_argument('-c', dest='config', type=str, help='config file') 
    subparser = parser.add_subparsers(title='subcommands', dest='subcommand')
    
    # parse commands related to design
    design_parser = subparser.add_parser('design', help='suboptions for gRNA design')
    design_parser.add_argument('-material', dest='material', default='rna', type=str, help='model material to use')
    design_parser.add_argument('-T', dest='temperature', default=37, type=float, help='temperature in C')
    design_parser.add_argument('-fstop', dest='fstop', default=0.05, type=float, help='minimum target defect')
    design_parser.add_argument('-s', dest='spec', type=str, help='cgRNA design specification to use')
    design_parser.add_argument('-d', dest='dim', nargs='+', type=int, help='domain dimensions')
    design_parser.add_argument('-gin', dest='gin', default='', type=str, help='csv file of RNA or protein sequences to use as trigger RNAs')
    design_parser.add_argument('-gout', dest='gout', default='', type=str, help='csv file of DNA sequences to use as guide RNA targets')
    design_parser.add_argument('-N', dest='Northo', default=1, type=int, help='number of orthogonal cgRNAs to design')
    design_parser.add_argument('-o', dest='outname', type=str, help='output file header')
    design_parser.add_argument('-scan', dest='scan', nargs=2, type=int, help='generate cgRNA library that scans an input RNA sequence with [window_size, stride]')
    
    # parse commands related to restarting design from checkpoint
    checkpoint_parser = subparser.add_parser('checkpoint', help='suboptions to restart designs from checkpoints')
    checkpoint_parser.add_argument('-trials', dest='trials', default=1, type=int, help='number of parallel designs to run')
    checkpoint_parser.add_argument('-cint', dest='ckpt_int', default=15, type=int, help='checkpoint interval in seconds')
    checkpoint_parser.add_argument('-maxopt', dest='max_opt', default=172800 ,type=float, help='max time in seconds to run the design')
    checkpoint_parser.add_argument('-i', dest='infile', nargs='+', type=str, help='nupack design file to load')
    checkpoint_parser.add_argument('-bucket', dest='bucket', type=str, help='s3 folder bucket for file upload and download')
    checkpoint_parser.add_argument('-url', dest='endpoint_url', type=str, help='s3 url')
    checkpoint_parser.add_argument('-id', dest='aws_access_key_id', type=str, help='s3 user id')
    checkpoint_parser.add_argument('-key', dest='aws_secret_access_key', type=str, help='s3 key')
    
    # parse commands related to file transfer
    ftp_parser = subparser.add_parser('ftp', help='suboptions for interacting with s3 fileserver')
    ftp_parser.add_argument('-push', dest='push', type=str, nargs='+', help='push files to s3 fileserver')
    ftp_parser.add_argument('-get', dest='get', type=str, nargs='+', help='get files containing keyword from s3 fileserver')
    ftp_parser.add_argument('-rm', dest='delete', type=str, nargs='+', help='delete files containing keyword from s3 fileserver')
    ftp_parser.add_argument('-ls', dest='ls', type=str, nargs='+', help='list files on the s3 fileserver')

    # parse commands related to k8 cluster
    kube_parser = subparser.add_parser('kube', help='suboptions for interacting with k8s')
    kube_parser.add_argument('-o', dest='outname', type=str, help='output file')
    kube_parser.add_argument('-clear', dest='clear', action='store_true', help='clear completed jobs')
    kube_parser.add_argument('-ls', dest='ls', default='*', type=str, help='list jobs')
    kube_parser.add_argument('-rm', dest='rm', type=str, help='remove jobs')
    kube_parser.add_argument('-i', dest='infile', nargs='+', type=str, help='input list of jobs to submit')
    kube_parser.add_argument('-pods', dest='pods', action='store_true', help='show pods')
    
    # parse commands related to analysis
    analysis_parser = subparser.add_parser('analysis', help='suboptions for strand analysis')
    analysis_parser.add_argument('-i', dest='infile', type=str, nargs='+', help='input files')
    analysis_parser.add_argument('-o', dest='outname', type=str, help='output file')
    analysis_parser.add_argument('-m', dest='method', type=str, help='analysis method to perform')
    
    analysis_parser.add_argument('-T', dest='temperature', type=float, default=37, help='temperature in C')
    analysis_parser.add_argument('-material', dest='material', type=str, default='rna', help='type of material')
    analysis_parser.add_argument('-s', dest='strands', nargs='+', help='strands to examine')
    analysis_parser.add_argument('-r', dest='reference', nargs='+', help='reference sequences to pass')
    analysis_parser.add_argument('-d', dest='dim', nargs='+', type=float, help='dimensions to pass')
    analysis_parser.add_argument('-mask', dest='mask', nargs='+', help='remove sequences from strands')
    analysis_parser.add_argument('-noterm', dest='noterm', action='store_true', help='remove terminator sequences from strands')
    analysis_parser.add_argument('-ex', dest='excludes', nargs='+', help='sequences to exclude from strands')
    analysis_parser.add_argument('-maxsize', dest='maxsize', type=int, help='max size of complexes')
    
    # parse commands related to unpacking results
    unpack_parser = subparser.add_parser('unpack', help='suboptions for unpacking nupack results')
    unpack_parser.add_argument('-i', dest='infile', type=str, nargs='+', help='input files')
    unpack_parser.add_argument('-o', dest='outname', type=str, help='output prefix')
    unpack_parser.add_argument('-d', dest='data', type=str, default='strands', choices=['strands','domains','tubes'], help='data type to extract')
    
    # parse commands related to oligos
    oligos_parser = subparser.add_parser('oligos', help='suboptions for oligos')
    oligos_parser.add_argument('-i', dest='infile', nargs='+', type=str, help='input csv files with strands to process')
    oligos_parser.add_argument('-o', dest='outfile', type=str, help='output filename for primers or oligos')
    oligos_parser.add_argument('-s', dest='strands', nargs='+', type=str, help='header of strand name to process')
    oligos_parser.add_argument('-m',dest='method', type=str, default='twist_outer', help='oligo synthesis method')
    oligos_parser.add_argument('-g', dest='ggsite', nargs='+', type=str, help='oligos overhang')
    oligos_parser.add_argument('-pad', dest='padding', type=int, default=300, help='pad ends with random base pairs so total length is at least N base pairs')
    oligos_parser.add_argument('-noterm', dest='noterm' , action='store_true', help='strip terminator sequences')
    oligos_parser.add_argument('-mask', dest='mask' , type=str, nargs='+', help='strip custom sequences from the oligos')
    oligos_parser.add_argument('-nonum', dest='renumber', action='store_false', help='do not renumber oligos')

    # parse arguments
    args = parser.parse_args()
    # initialize the log file
    set_log()
    
    if args.verbose:
        print(args)
    
    # options to generate nupack design specs
    if args.subcommand=='design':
        opt = Designer(args)
        opt.run()

    # options for starting jobs from design specs
    elif args.subcommand=='checkpoint':
        opt = Worker(args)
        opt.run()

    # options for interfacing with s3 fileserver
    elif args.subcommand=='ftp':
        server = Fileserver(args)
        if args.ls!=None:
            files = server.get_filelist(args.ls)
            logging.info(str(np.array(files)))
        elif args.get!=None:
            files = server.get_filelist(args.get)
            server.fetch_files(files)
        elif args.push!=None:
            server.upload_files(args.push)
        elif args.delete!=None:
            files = server.get_filelist(args.delete)
            server.delete_files(files)
    
    # options for interfacing with kubernetes cluster
    elif args.subcommand=='kube':
        kctrl = Kube(args)
        if args.infile!=None:
            kctrl.submit_jobs()
        elif args.pods:
            if args.clear:
                kctrl.clear_pods()
            elif args.rm!=None:
                kctrl.delete_pods()
            else:
                print(kctrl.list_pods())
        elif args.rm!=None:
            kctrl.delete_jobs()
        elif args.clear:
            kctrl.clear_completed() 
        else:
            print(kctrl.list_jobs())
            
    # options for unpacking nupack results
    elif args.subcommand=='unpack':
        parse_results(args.infile, args.outname, args.data)
        
    # run strand analysis
    elif args.subcommand=='analysis':
        opt = Analysis(args)
        opt.run()
        
    # options for generating oligos
    elif args.subcommand=='oligos':
        opt = Oligos(args)
        opt.run() 

if __name__ == '__main__':
    main()
