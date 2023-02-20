#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# Class function for manipulating the kubernetes cluster

import kubernetes as k8s
#import kubernetes.client
import pandas as pd
import numpy as np

from .fileserver import Fileserver
from .misc import get_timestamp
from .misc import get_hash
from .misc import update_parameters
from .misc import read_config

import hashlib
import datetime
import time
import sys
import fnmatch
import logging

class Kube:
    params = {}
    params['k8'] = {'host':'',
                    'api_key':'',
                    'api_key_prefix':'',
                    'namespace':'',
                    'container_name':'',
                    'container_image':'',
                    'regcred':'',
                    'cpu':['5000m','5000m'],
                    'memory':['16G','16G']}
    params['s3'] = {'bucket':'',
                    'endpoint_url':'',
                    'aws_access_key_id':'',
                    'aws_secret_access_key':''}

    def __init__(self, args=None):
        self.params = update_parameters(args, self.params)
        # read config if available
        if 'config' in self.params:
            data = read_config(self.params['config'])
            self.params = update_parameters(data, self.params)
            del(self.params['config'])
        # test k8 credentials
        self.set_credentials()
        self.test_credentials()
        self.server = Fileserver(self.params)
        
    def set_credentials(self):
        '''
        This function sets the BatchV1Api and CoreV1Api for the kube class
        '''
        params = self.params['k8']
        self.configuration = k8s.client.Configuration()
        self.configuration.api_key['authorization'] = params['api_key']
        self.configuration.api_key_prefix['authorization'] = params['api_key_prefix']
        self.configuration.host = params['host']
        if 'ssl_cert' in params.keys():
            self.configuration.ssl_ca_cert=params['ssl_cert']
        # apply config to get API
        self.BatchV1Api = k8s.client.BatchV1Api(k8s.client.ApiClient(self.configuration))
        self.CoreV1Api = k8s.client.CoreV1Api(k8s.client.ApiClient(self.configuration))

    def test_credentials(self):
        '''
        Function to test kubernetes credentials
        '''
        try:
            api_response = self.BatchV1Api.get_api_resources()
            if self.params['verbose']:
                print(api_response)
        except k8s.client.rest.ApiException as e:
            print('Exception when calling API: %s\n' % e)
            sys.exit(1)
    
    def list_pods(self):
        '''
        Lists the pods currently running on the cluster
        returns list of pods as a pandas dataframe
        '''
        try:
            ret = self.CoreV1Api.list_namespaced_pod(self.params['k8']['namespace'])
            out = []
            for i in ret.items:
                out.append([i.status.pod_ip, i.metadata.namespace, i.metadata.name, i.status.phase])
            out = pd.DataFrame(out, columns=['ip','namespace','name','phase'])
            # filter for pods of interest
            key = self.params['ls']
            jobs = fnmatch.filter(out['name'].values, key)
            out = out[out['name'].isin(jobs)]
            if 'outname' in self.params:
                out.to_csv(self.params['outname'], index=False)
            return out
        except k8s.client.rest.ApiException as e:
            print('Exception when calling API: %s\n' % e)
            sys.exit(1)

    def delete_pods(self, pods=[]):
        '''
        Deletes pods on kubernetes that match a list of names or keyword
        pods = list of pod names to delete
        '''
        if 'rm' in self.params:
            pods = self.list_pods()
            key = self.params['rm']            
            jobs = fnmatch.filter(pods['name'].values, key)
            pods = pods[pods['name'].isin(jobs)]
 
        for namespace, name, phase in pods[['namespace','name','phase']].values:
            try:
                print('deleting pod',namespace, name, phase)
                api_response = self.CoreV1Api.delete_namespaced_pod(name, namespace)
                if self.params['verbose']:
                    print(api_response)
            except k8s.client.rest.ApiException as e:
                print('Exception when calling API: %s\n' % e)
                sys.exit(1)

    def clear_pods(self):
        '''
        Clears pods that errored out or finished runnng
        '''
        pods = self.list_pods()
        pods = pods[~(pods['phase'] == 'Running')&~(pods['phase'] == 'Pending')]
        self.delete_pods(pods)

    def list_jobs(self):
        '''
        Returns list of currently running jobs on kubernetes
        '''
        key = '*'
        if 'ls' in self.params:
            key = self.params['ls']
        ret = self.BatchV1Api.list_namespaced_job(self.params['k8']['namespace'])
        out = []
        for x in ret.items:
            job_name = x.metadata.name
            start_time = x.status.start_time
            end_time = x.status.completion_time
            active = x.status.active
            failed = x.status.failed
            succeeded = x.status.succeeded
            cmd = ''
            out.append([job_name, start_time, end_time, failed, succeeded, cmd])
        col = ['name','start','end','failed','succeeded','cmd']
        out = pd.DataFrame(out, columns=col)
        jobs = fnmatch.filter(out['name'].values, key)
        out = out[out['name'].isin(jobs)]
        if 'outname' in self.params:
            out.to_csv(self.params['outname'], index=False)
        return out.sort_values(by=['name'])

    def delete_jobs(self):
        '''
        Deletes jobs matching a keyword
        '''
        namespace = self.params['k8']['namespace']
        key = self.params['rm']
        self.params['ls'] = key
        df = self.list_jobs()
        for name, start, end, failed, succeeded in df[['name','start','end','failed','succeeded']].values:
            print('deleting', name, start, end, failed, succeeded)
            ret = self.BatchV1Api.delete_namespaced_job(name=name, namespace=namespace)
            if self.params['verbose']:
                print(ret)

    def clear_completed(self):
        '''
        Deletes completed jobs from the queue
        '''
        namespace = self.params['k8']['namespace']
        df = self.list_jobs()
        df = df[df['end'].notna()]
        for name, start, end, failed, succeeded in df[['name','start','end','failed','succeeded']].values:
            print('deleting', name, start, end, failed, succeeded)
            ret = self.BatchV1Api.delete_namespaced_job(name=name, namespace=namespace)
            if self.params['verbose']:
                print(ret)

    def create_job(self, name, cmd, env_vars={}):
        '''
        Adds a job onto the kubernetes cluster
        '''
        params = self.params['k8']
        namespace = params['namespace']
        container_image = params['container_image']
        container_name = params['container_name']
        regcred = params['regcred']

        # Body is the object Bodyv
        body = k8s.client.V1Job(api_version='batch/v1', kind='Job')
        # Body needs Metadata
        # Attention: Each JOB must have a different name!
        body.metadata = k8s.client.V1ObjectMeta(namespace=namespace, name=name)
        # And a Status
        body.status = k8s.client.V1JobStatus()
        # Now we start with the Template...
        template = k8s.client.V1PodTemplate()
        template.template = k8s.client.V1PodTemplateSpec()
        # set cpu and memory requirements
        req = {}
        limit = {}
        for key in ['cpu','memory']:
            req[key] = params[key][0]
            limit[key] = params[key][1]
        resources = k8s.client.V1ResourceRequirements(limits=limit, requests=req)
        # add any environment variables
        env_list = []
        for env_name, env_value in env_vars.items():
            env_list.append(k8s.client.V1EnvVar(name=env_name, value=env_value))
        # specify the container
        container = k8s.client.V1Container(name=container_name, image=container_image, env=env_list, command=cmd, resources=resources)
        #toleration = k8s.client.V1Toleration(effect='NoSchedule', key='nupack.org/nupack-worker', operator='Exists', value='nupack-prod')
        toleration = k8s.client.V1Toleration(effect='NoSchedule', key='nupack.org/nupack-worker', operator='Exists')
        template.template.spec = k8s.client.V1PodSpec(containers=[container], tolerations=[toleration], restart_policy='Never', image_pull_secrets=[k8s.client.V1LocalObjectReference(regcred)])

        # And finaly we can create our V1JobSpec!
        body.spec = k8s.client.V1JobSpec(template=template.template)
        try: 
            api_response = self.BatchV1Api.create_namespaced_job(namespace, body, pretty=True)
            if self.params['verbose']:
                print(api_response)
        except k8s.client.rest.ApiException as e:
            print('Exception when calling API->create_namespaced_job: %s\n' % e)

    def submit_jobs(self):
        '''
        Submits jobs to the kubernetes cluster
        '''
        files = self.params['infile']
        files = np.array(files)
        c = np.array([('.spec' in f) for f in files])
        spec = [f for f in files[c]]
        bash = files[~c]
        for f in bash:
            try:
                self.submit_bash(f)
            except:
                print('Failed to submit bash file:',f)
        if len(spec) > 0:
            self.submit_nupack(spec)

    def submit_bash(self, infile):
        '''
        Submit a bash script
        '''
        # upload bash script
        self.server.upload_files(infile)
        # make the name
        jname = Kube.format_name(infile)
        # submit the job
        print('submitting '+jname)
        cmd = 'CGDesigner checkpoint'
        cmd+= ' -url '+str(self.params['s3']['endpoint_url'])
        cmd+= ' -id '+str(self.params['s3']['aws_access_key_id'])
        cmd+= ' -key '+str(self.params['s3']['aws_secret_access_key'])
        cmd+= ' -bucket '+str(self.params['s3']['bucket'])
        cmd+= ' -i '+infile
        print(cmd)
        self.create_job(jname, cmd.split(' '))

    def submit_nupack(self, files):
        '''
        Submit nupack jobs
        '''
        # upload files to s3
        self.server.upload_files(files)
        trials = self.params['trials']
        cmd = 'CGDesigner checkpoint'
        cmd+= ' -cint '+str(self.params['ckpt_int'])
        cmd+= ' -maxopt '+str(self.params['max_opt'])
        cmd+= ' -trials '+str(self.params['trials'])        
        cmd+= ' -url '+str(self.params['s3']['endpoint_url'])
        cmd+= ' -id '+str(self.params['s3']['aws_access_key_id'])
        cmd+= ' -key '+str(self.params['s3']['aws_secret_access_key'])
        cmd+= ' -bucket '+str(self.params['s3']['bucket'])
        cmdbase = cmd

        files = np.array_split(files, 1+int(len(files)/trials))
        for batch in files:
            name = batch[0].split('.spec')[0].split('/')[-1]
            name = self.params['s3']['bucket'] + '_' + name
            jname = Kube.format_name(name)
            if jname[-1]=='-':
                jname = jname[:-1]
            # submit the job
            print('submitting '+jname)
            cmd = cmdbase + ' -i '+' '.join(batch)
            print(cmd)
            self.create_job(jname, cmd.split(' '))

    def format_name(name):
        h = get_hash(name)
        jname = h[:10]+'_'+name
        jname = jname.replace('_','-')
        jname = jname.replace(' ','')
        jname = jname.replace('.','-')
        jname = jname.lower()
        jname = jname[:63]
        return jname


