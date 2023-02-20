#!/usr/bin/env python
# Generate oligos for cloning

import Bio
import pandas as pd
import numpy as np
import fnmatch
import logging

from .formula import Formula
from .designer import Designer
from .misc import update_parameters

class Oligos:
    params = {}
    def __init__(self, args=None):
        # defines extraneous sequences to remove from analysis like terminators
        parts = Formula().params['parts']
        self.params['mask'] = [parts[k]['sequence'].lower() for k in parts.keys()]
        # defines sequences to exclude
        parts = Designer().params
        x = [parts['enzymes'][k] for k in parts['enzymes'].keys()]
        self.params['excludes'] = x
        self.params['excludes']+= [str(Bio.Seq.Seq(i).reverse_complement()) for i in x]
        self.params['excludes']+= parts['excludes']
        self.params['excludes'] = [k.lower() for k in self.params['excludes']]        
        # load user defined parameters
        self.params = update_parameters(args, self.params)

    def run(self):
        '''
        Get oligos for cloning
        '''
        out = []
        for fname in self.params['infile']:
            print('generating oligos for '+fname)
            df = pd.read_csv(fname)
            df = Oligos.get_tagged(df, self.params['strands'])
            # apply proper descriptors
            if 'filename' in df.columns:
                name = [f.split('/')[-1].split('_trial')[0].split('_d0.')[0] for f in df.filename]
            else:
                name = fname.split('/')[-1].split('_trial')[0].split('_d0.')[0]
            df['filename'] = name
            if ('description' in df.columns) == False:
                df['description'] = df['filename']+'_'+df['Strand']
            out.append(df)
        df = pd.concat(out)    

        # strip masked sequences
        df['Sequence'] = Oligos.strip_sequence(df.Sequence.values, self.params['mask'])
    
        # run oligo synthesis method
        tube_func = getattr(self, str(self.params['method']))
        df = tube_func(df)
        # make all columns lower case
        for c in df.columns:
            df = df.rename(columns={c:c.lower()})
        
        if self.params['renumber']:
            df['Strand'] = [str(i)+'_'+df.iloc[i]['description'] for i in range(len(df))]
        else:
            df['Strand'] = df['description']

        df['Strand'] = [i[:30] for i in df['Strand']]
        df['length'] = [len(i.replace(' ','')) for i in df['sequence']]
        df = df.rename(columns={'sequence':'Sequence'})
        col = ['length','description','Strand','Sequence',]
        df[col].to_csv(self.params['outfile'],index=False)
        return df

    def get_tagged(df, tags):
        '''
        Keep only entries matching tags
        '''
        m = []
        for k in tags:
            m+= list(fnmatch.filter(df['Strand'].values, k))
        m = list(set(m))
        return df[df['Strand'].isin(m)]

    def twist_outer(self, df):
        '''
        Pad outside only
        '''
        seqlist = self.twist(df)
        out = []
        for seq in seqlist:
            s = self.add_padding([''.join(seq[1:])], self.params['padding'])
            out.append([seq[0], s])
        col = ['description','Sequence']
        return pd.DataFrame(out, columns=col)

    def twist_inner(self, df):
        '''
        Pad between sequences
        '''
        seqlist = self.twist(df)
        out = []
        for seq in seqlist:
            s = self.add_padding(seq[1:], self.params['padding'])
            out.append([seq[0], s])
        col = ['description','Sequence']
        return pd.DataFrame(out, columns=col)
 
    def twist(self, df):
        # group strands based on orthogonal sets
        seqlist = Oligos.get_orthoset(df, self.params['strands'])
        # get ggsites
        if '.csv' in self.params['ggsite'][0]:
            ggsite = pd.read_csv(self.params['ggsite'][0])
            ggsite = ggsite[ggsite.columns[ggsite.columns!='label']]
            ggsite = ggsite.values
        else:
            ggsite = np.array([self.params['ggsite']])
        
        # add ggsites
        N = len(ggsite)
        for j in range(1, len(seqlist[0])):
            for i in range(len(seqlist)):
                seqlist[i][j] = ggsite[i%N,(j-1)*2] + seqlist[i][j] + ggsite[i%N,(j-1)*2+1]
        return seqlist

    def twist_pool(self, df):
        '''
        Generate oligos
        '''
        # get ggsites
        if '.csv' in self.params['ggsite'][0]:
            dfg = pd.read_csv(self.params['ggsite'][0])
 
        dims = self.params['offset']
        if sum(dims) > len(dfg):
            logging.error('Not enough primers sites. dims='+str(dims)+' sites='+str(len(dfg)))
            sys.exit(1)

        # group strands based on orthogonal sets
        df['idx'] = [i for i in range(len(df))]
        for i in range(len(dims)):
            d = dims[i]
            # primers to use
            x = df['idx']%d
            x = x.astype(int)
            df['idx'] = (df['idx'] - x)/d
            # add primers
            pmr = dfg.seqL + dfg.primer + dfg.seqR
            pmr = pmr.values
            if i%2==0:
                df['Sequence'] = [pmr[j]+' ' for j in x] + df['Sequence']
            else:
                df['Sequence'] = df['Sequence'] + [' '+pmr[j] for j in x]
            # add names 
            label = dfg.label.astype(str).values
            df['description'] = [label[j]+'_' for j in x] + df['description']
            # use next set of primers
            dfg = dfg.iloc[d:]
        return df

    def strip_sequence(seq, mask):
        '''
        Strip sequences from oligos
        '''
        # lower case everything
        mask = [k.lower().replace('u','t') for k in mask]
        seq = [s.lower() for s in seq]
        # strip sequences 
        for k in mask:
            seq = [s.replace(k,'') for s in seq]
        return seq

    def get_orthoset(df, strands):
        '''
        Group strands into orthogonal sets
        '''
        files = df['filename'].drop_duplicates().values
        filename = df['filename'].values
        name = df['Strand'].values
        seq = df['Sequence'].values

        out = []
        for f in files:
            fname = name[filename==f]
            fseq = seq[filename==f]
            tag = np.array([i.split('[')[1].split(']')[0] for i in fname])
            utag = np.unique(tag)
            for k in utag:
                sname = fname[tag==k]
                slist = fseq[tag==k]
                outname = f
                outseq = []
                for s in range(len(strands)):
                    m = fnmatch.filter(sname, strands[s])
                    loc = [(m[0] in n) for n in sname]
                    outname+= m[0]+'_'
                    outseq.append(slist[loc][0])
                out.append([outname[:-1]] + outseq)
        return out
    
    def get_rand_seq(self, length, alphabet='ATGC'):
        '''
        Generate random sequence
        '''
        idx = np.random.randint(0, high=len(alphabet), size=length)
        return ''.join([alphabet[i] for i in idx])

    def add_padding(self, seqlist, length, maxiter=1000):
        '''
        Find new padding sequence until exclude sequences are not found
        '''
        exclude = [i for i in list(self.params['excludes'])]
        # if total length is good, then 
        out = ''.join(seqlist)
        if len(out) >= length:
            return out
        # otherwise pad the space between sequences
        N = int((length - len(out))/(len(seqlist)+1)) + 1
        for k in range(maxiter):
            out = ''.join([self.get_rand_seq(N) + seqlist[i] for i in range(len(seqlist))])
            out+= self.get_rand_seq(N)
            # check for gg sites
            ggsite = 0
            for j in exclude:
                ggsite+=self.check_sequence(j, out)
                if ggsite > 0:
                    break
            if ggsite==0:
                break
        return out
 
    def check_sequence(self, query, target, rev_comp=True):
        '''
        Check if query sequence in contained within target sequence
        '''
        queryR = str(Bio.Seq.Seq(query).reverse_complement())
        queryR = queryR.lower()
        target = target.lower()
        return (query in target) + (queryR in target)
