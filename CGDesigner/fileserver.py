#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# Class function for interfacing with the s3 fileserver

import numpy as np
import botocore
import boto3

import fnmatch
import logging
from .misc import update_parameters
from .misc import read_config

class Fileserver:
    params = {}
    def __init__(self, args=None):
        self.params = update_parameters(args, self.params)
        # read config if available
        if 'config' in self.params:
            data = read_config(self.params['config'])
            self.params = update_parameters(data, self.params)
            del(self.params['config'])
        self.client = self.get_client()
        self.server = self.get_fileserver()

    def get_fileserver(self):
        '''
        Configure s3 file server
        '''
        try:
            return boto3.resource('s3',
                    endpoint_url=self.params['s3']['endpoint_url'],
                    aws_access_key_id=self.params['s3']['aws_access_key_id'],
                    aws_secret_access_key=self.params['s3']['aws_secret_access_key'],
                    config=botocore.client.Config(signature_version='s3v4'),
                    region_name='us-east-1')
        except:
            logging.info('Failed to initialize fileserver')
            return None

    def get_client(self):
        '''
        Configure s3 file server
        '''
        try:
            return boto3.client(service_name='s3',
                    endpoint_url=self.params['s3']['endpoint_url'],
                    aws_access_key_id=self.params['s3']['aws_access_key_id'],
                    aws_secret_access_key=self.params['s3']['aws_secret_access_key'],
                    config=botocore.client.Config(signature_version='s3v4'),
                    region_name='us-east-1')
        except:
            logging.info('Failed to initialize fileserver client')

    def get_filelist(self, keys=['']):
        '''
        Get list of files in s3 bucket
        '''
        try:
            print('Searching for ', keys)
            server = self.client

            files = []
            for key in keys:
                prefix = ''
                if '/' in key:
                    prefix = key.split('/')[0]
                elif key[0]!='*':
                    prefix = key.split('*')[0]       
                files+= [i['Key'] for i in server.list_objects(Bucket=self.params['s3']['bucket'], Prefix=prefix)['Contents']]
            files = [f for f in np.unique(files)]

            print('total files', len(files))
            keep = []
            for k in keys:
                keep += fnmatch.filter(files, k)
            keep = np.array(list(set(keep)))
            print('files matching keys ', len(keep))
            # remove duplicates
            s = np.argsort(keep)
            return keep[s]
        except:
            logging.info('Failed to search for files on fileserver')
            return []

    def fetch_files(self, files):
        '''
        Fetch files from s3 fileserver
        '''
        server = self.server
        for fname in files:
            try:
                print('Downloading '+fname)
                server.Bucket(self.params['s3']['bucket']).download_file(Key=fname, Filename=fname)
            except:
                logging.info('Failed to fetch '+fname)

    def upload_files(self, files):
        '''
        Upload files to s3 fileserver 
        '''
        server = self.server
        if type(files) == str or type(files)==np.str_:
            files = [files]
        for f in files:
            print('uploading '+f)
            infile = f
            outfile = f
            try:
                server.meta.client.upload_file(Filename=infile, Bucket=self.params['s3']['bucket'], Key=outfile)
            except:
                logging.info('Failed to upload '+infile)

    def delete_files(self, files):
        '''
        Delete files from s3 fileserver
        '''
        server = self.server
        for fname in files:
            print('Deleting '+fname+' from s3 fileserver bucket '+self.params['s3']['bucket'])
            try:
                server.Object(self.params['s3']['bucket'], fname).delete()
            except:
                logging.info('Failed to delete '+fname)



