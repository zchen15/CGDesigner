#!/usr/bin/env python
# analysis of RNA/DNA strands

# computing libraries
import nupack
import numpy as np
import scipy as sp
import pandas as pd
import Bio
import Bio.Seq

# plotting libraries
import matplotlib.pyplot as plt
import logging

# custom libraries
from .formula import Formula
from .designer import Designer
from .oligos import Oligos
from .misc import update_parameters
from .misc import read_csv

class Analysis:
    params = {'kB':0.001985875}

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
        self.params['beta'] = 1/(self.params['kB']*(self.params['temperature']+273.15))
        self.model = nupack.Model(material=self.params['material'], celsius=self.params['temperature'])

    def run(self):
        '''
        Selects which function in the obj to run
        '''
        if self.params['verbose']:
            logging.info(str(self.params))
        df = self.load_data()
        tube_func = getattr(self, str(self.params['method']))
        df = tube_func(df)
        df.to_csv(self.params['outname'], index=False)

    def load_data(self):
        '''
        Standardized data loading for this object
        only get tagged strands
        strip terminator sequences
        '''
        df = read_csv(self.params['infile'])
        # get only tagged sequences
        if 'strands' in self.params.keys():
            df = Oligos.get_tagged(df, self.params['strands'])
        # strip terminator sequence
        if self.params['mask']!=[]:
            df['Sequence'] = Oligos.strip_sequence(df['Sequence'].values, self.params['mask'])
        return df

    def remove_toeholds(self, df):
        '''
        Extract toehold sequences of toehold switch
        '''
        # extract first 30nt
        df['Sequence'] = [i[:30] for i in df['Toehold']]
        # align against reference
        df= self.match_to_reference(df)
        return df

    def remove_homopolymers(self, df):
        '''
        Remove homopolymers from results list
        '''
        # exclusion list
        excludes = self.params['excludes']
        # mark homopolymers
        name = df['filename']+'_'+df['Strand']
        c = Analysis.get_excluded(df['Sequence'].values, excludes, name=name.values)            
        keep = df[~c]['filename'].drop_duplicates().values
        # reload data
        del(self.params['strands'])
        df = self.load_data()
        # keep strands without homopolymers
        return df[df['filename'].isin(keep)]
    
    def get_excluded(seq, ex, name):
        '''
        Returns numpy array of true/false if entries are part of exclusion list
        '''
        toss = np.array([False]*len(seq))
        for i in range(len(seq)):
            for x in ex:
                c = (x.lower() in seq[i].lower())
                toss[i]+= c
                if c:
                    logging.info(x+' found in '+name[i])
        return toss

    def match_to_reference(self, df):
        '''
        Match a sequence to a reference strand
        Return dataframe with reference, pos1, pos2, orientation, length columns
        '''
        ref = read_csv(self.params['reference'])
        out = []
        for refn, refseq in ref[['Strand','Sequence']].values:
            logging.info('aligning strands against '+refn)
            loc1 = []
            loc2 = []
            for seq in df['Sequence']:
                # search in fwd
                x = refseq.lower().find(seq.lower())
                loc1.append(x)
                # search in reverse
                revseq = str(Bio.Seq.Seq(seq).reverse_complement())
                x = refseq.lower().find(revseq.lower())
                loc2.append(x)
            df.loc[:,'reference'] = refn
            df1 = df.copy()
            df2 = df.copy()
            df1['pos1'] = loc1
            df1['orientation'] = 1
            df2['pos1'] = loc2
            df2['orientation'] = -1
            out.append(df1)
            out.append(df2)
        out = pd.concat(out)
        out = out[out['pos1'] > -1]
        out['length'] = [len(i) for i in out['Sequence']]
        out['pos2'] = out['pos1']+out['length']
        return out
        
    def get_ref_window(self, df):
        '''
        Find location on reference and add window of sequence around it to dataframe
        Returns dataframe with extra row for each trigger that has N base pairs on each side
        '''
        df = self.match_to_reference(df)
        width = int(self.params['dim'][0])
        ref = read_csv(self.params['reference'])
        out = []
        for refn, refseq in ref[['Strand','Sequence']].values:
            logging.info('assigning windows against '+refn)
            wseq = []
            c = df['reference'] == refn
            x = df[c]
            s0 = x['pos1'].values 
            s1 = s0 - width 
            s2 = s0 + x['length'].values + width
            s1[s1 < 0] = 0
            s2[s2 > len(refseq)] = len(refseq)
            x.loc[:,'Sequence'] = [refseq[s1[i]:s2[i]] for i in range(len(s1))]
            out.append(x)
        return pd.concat(out) 

    def get_dG(self, seq):
        '''
        Get free energy of a set of complexes
        '''
        strands = [nupack.Strand(seq[i], name='S'+str(i)) for i in range(len(seq))]
        c1 = nupack.Complex(strands, name='c1')
        set1 = nupack.ComplexSet(strands=strands, complexes=[c1])
        # calculate the partition function for each complex in the complex set
        res = nupack.complex_analysis(complexes=set1, model=self.model, compute=['pfunc'])
        dG = res[c1].free_energy
        return dG 

    def get_trigger_dG(self, df):
        '''
        Analyze binding energy of trigger strands
        '''
        out = []
        for seq in df['Sequence'].values:
            # single stranded
            dGss = nupack.pfunc(strands=[seq], model=self.model)
            # unbase paired single stranded
            dGup = nupack.structure_energy(strands=[seq], structure='.'*len(seq), model=self.model)
            # double stranded
            revseq = str(Bio.Seq.Seq(seq).reverse_complement())
            dGds = nupack.pfunc(strands=[seq, revseq], model=self.model)
            # full duplex
            dGdss = nupack.structure_energy(strands=[seq, revseq], structure='('*len(seq)+'+'+')'*len(seq), model=self.model)
            x = [dGss[1], dGup, dGds[1], dGdss]
            out.append(x)
        out = np.array(out)
        col = ['dGss','dGssu','dGds','dGdsu']
        for i in range(len(col)):
            df[col[i]] = out[:,i]
        return df

    def get_strand_window(df, length, stride):
        # get window size and stride
        length = int(length)
        stride = int(stride)
        # generate a list of strands
        strands = []
        for name, seq in df[['Strand','Sequence']].values:
            w1 = np.arange(0, len(seq), stride)
            w2 = w1 + length
            w2[w2 > len(seq)] = len(seq)
            strands+= [[name, w1[i], w2[i], seq[w1[i]:w2[i]]] for i in range(len(w1))]
        col = ['Strand','pos1','pos2','Sequence']
        out = pd.DataFrame(strands, columns=col)
        return out

    def get_mRNA_dG(self, df):
        '''
        get dG of mRNA windows
        '''
        length, stride = self.params['dim']
        df = Analysis.get_strand_window(df, length, stride)
        return self.get_trigger_dG(df)
 
    def get_mRNA_dGx(self, df):
        '''
        get cross bind dG of mRNA windows
        '''
        length, stride = self.params['dim']
        df = Analysis.get_strand_window(df, length, stride) 
        out = []
        c = 0
        seqlist = np.array()
        for g in np.unique(df['name'].values):
            x = df[df['name']==g]
            w = x['pos1']+'_'+x['pos2']
            w = w.values
            seq = x['sequence'].values 
            for i in range(len(w)):
                for j in range(i, len(w)):
                    dGx = nupack.pfunc(strands=[seq[i], seq[j]], model=self.model)
                    out.append([g, w[i], w[j], dGx[1], dGx[0]])
                    # track progress
                    if c%100==0:
                        logging.info('calculating '+g+' for '+str(i)+','+str(j)+' of '+str(len(w))+' strands')
                    c+=1
        return pd.DataFrame(out, columns=['name','w1','w2','dGx','Q'])

    def parse_defect(self, df):
        '''
        Parse design defect from the filename
        '''
        df['defect'] = [float(i.split('_d0')[-1].split('_strand')[0]) for i in df['filename']]
        return df

    def add_position(self, df):
        '''
        Add trigger position information to the dataframe
        '''
        trig = self.match_to_reference(df)
        # reload data
        del(self.params['strands'])
        df = self.load_data()
        
        col = ['filename','reference','pos1','pos2','length','orientation']
        for k in col[1:]:
            if k in df.columns:
                df = df.drop(columns=[k])
        df = df.merge(trig[col], on='filename', how='left')
        df = df.dropna()
        return df
    
    def add_prediction(self, df):
        '''
        Add prediction score about trigger position
        '''
        # add prediction score
        ref = read_csv(self.params['reference'])
        pcol = 'prediction'
        out = []
        for t in ref['reference'].drop_duplicates().values:
            dfr = ref[ref['reference'] == t]
            dfd = df[df['reference'] == t]
            if len(dfr) > 0 and len(dfd) > 0:    
                x1 = dfr['pos1'].values
                x2 = dfr['pos2'].values
                y1 = dfr[pcol].values
                x,y = Analysis.get_window_average(x1, x2, y1)
                pred = []
                for pos1, pos2 in dfd[['pos1','pos2']].values:
                    idx = (x >= pos1) & (x < pos2)
                    pred.append(np.mean(y[idx]))
                dfd[pcol] = pred
                out.append(dfd)
        return pd.concat(out)

    def get_n_best(self, df):
        '''
        Sort to N best non-overlapping designs using cost function
        '''
        nsort, L = self.params['dim']
        nsort = int(nsort)
        # auto generate prediction
        mcol = self.params['reference'][0]
        if mcol == 'defect':
            df['defect'] = 1/df['defect']
        
        # get sort for the best designs
        cols = ['filename','pos1','pos2',mcol]
        cols = np.unique(cols)
        x = df[cols].drop_duplicates()
        s = np.argsort(x[mcol].values)[::-1]
        
        # find best non-overlapping designs
        y = x.iloc[s]
        s = Analysis.get_overlappers(y['pos1'].values, y['pos2'].values, L=L)
        y = y.iloc[~s]
        y['rank'] = [i for i in range(len(y))]
        
        # merge ranking information
        if 'rank' in df.columns:
            df = df.drop(columns=['rank'])
        col = ['filename','rank']
        df = df.merge(y[col], on='filename', how='left')
        
        # set a cut off
        total = len(df['filename'].drop_duplicates())
        df = df[df['rank'] < nsort]
        N = len(df['filename'].drop_duplicates())
        logging.info('filtered '+str(N)+'/'+str(total))
        
        # reverse defect computation
        if mcol == 'defect':
            df['defect'] = 1/df['defect']
        return df

    def get_overlappers(pos1, pos2, L=0):
        '''
        Keep only first set of non-overlapping results in list
        '''
        N = len(pos1)
        keep = np.array([True]*N)
        for i in range(0,N):
            if keep[i]:
                for j in range(i+1,N):
                    if keep[j]:
                        x1 = [pos1[i], pos2[i]]
                        x2 = [pos1[j], pos2[j]]
                        keep[j] = ~Analysis.is_overlapping(x1, x2, L)
        return ~keep 

    def is_overlapping(x1, x2, L=0):
        '''
        Check if pairs of coordinates are overlapping
        '''
        # x1 > x2
        if x1[1] > x2[1]:
            return Analysis.is_overlapping(x2, x1, L)
        elif (L > 0 and L < 1):
            # flag based on threshold overlap
            d = x1[1] - x2[0]
            L1 = x1[1] - x1[0]
            L2 = x2[1] - x2[0]
            return L < d/min(L1, L2)
        else:
            # flag overlap based on distance
            return x2[0] - x1[1] < L

    def on_target_analysis(self, df):
        '''
        Run on target tube analysis
        '''
        # load strand info
        conc = np.array(self.params['dim'])*1e-9
        n_strands = len(np.unique(df['Strand']))
        conc = [conc[i] for i in range(n_strands)]
        if ('Concentration (M)' in df.columns)==False:
            df['Concentration (M)'] = [i for i in conc]*int(len(df)/len(conc))

        spec = nupack.SetSpec(max_size=self.params['maxsize'])
        out = []
        for fname in np.unique(df['filename']):
            logging.info('analyzing '+fname)
            c = df['filename']==fname
            strands = [nupack.Strand(s, name=n) for n,s in df[c][['Strand','Sequence']].values]
            conc = [cnc for cnc in df[c]['Concentration (M)'].values]
            ts = {strands[i]:conc[i] for i in range(len(strands))}
            cp_set = nupack.ComplexSet(strands=strands, complexes=spec)
            # run complex analysis
            cp_res = nupack.complex_analysis(complexes=cp_set, model=self.model, compute=['pfunc'])
            # run tube analysis
            cc_res = nupack.complex_concentrations(tube=cp_set, data=cp_res, concentrations=ts)
            x = Analysis.get_dataframe(cp_res, cc_res)
            x['filename'] = fname
            out.append(x)
        return pd.concat(out)

    def get_dataframe(complex_result=None, complex_concentration=None):
        '''
        Get complex result and/or complex concentration as dataframe
        '''
        # convert to dataframe
        if type(complex_result) != type(None):
            complex_result = complex_result.complex_frame(pretty=False)
        if type(complex_concentration) != type(None):
            complex_concentration = complex_concentration.tube_frame(pretty=False)
            complex_concentration = complex_concentration.rename(columns={'tube':'Concentration (M)'})
        
        # merge dataframes
        if type(complex_result) != type(None) and type(complex_concentration) != type(None):
            df = complex_result.merge(complex_concentration, on='complex', how='left')
        elif type(complex_result) != type(None):
            df = complex_result
        elif type(complex_concentration) != type(None):
            df = complex_concentration
        else:
            return None
        # format complex to string
        df['complex'] = [i.name.replace('(','').replace(')','') for i in df['complex']]
        df = df.rename(columns={'complex':'Complex'})
        return df

    def off_target_analysis(self, df):
        '''
        Run off target tube analysis against mRNA window
        '''
        ref = read_csv(self.params['reference'])
        A, B, s1, L2, s2 = self.params['dim']
        [A, B, s1, L2, s2] = [int(i) for i in [A, B, s1, L2, s2]]
        # get ref sequence windows
        ref = Analysis.get_strand_window(ref, L2, s2)
        ref['name'] = ref['Strand'] + '_w' + ref['pos1'].astype(str) + '_' + ref['pos2'].astype(str)
        #ref['name'] = ref['Strand']
        # get query sequence windows
        if A+abs(B) > 0:
            query = Analysis.get_strand_window(df, A+abs(B), s1)
            # get reverse complement
            query['Sequence'] = [str(Bio.Seq.Seq(i).reverse_complement()) for i in query['Sequence']]
            #query['name'] = query['Strand']
            query['name'] = query['Strand'] + '_w' + query['pos1'].astype(str) + '_' + query['pos2'].astype(str)
            query['Strand'] = 'seq1'
            # remove sequences too short
            query['length'] = [len(i) for i in query['Sequence']]
            query = query[query['length'] == A+abs(B)]
            # generate hairpin
            '''
            if B > 0:
                query['Sequence'] = [str(Bio.Seq.Seq(i[:B]).reverse_complement()) + 'TTT' + i for i in query['Sequence']]
            elif B < 0:
                query['Sequence'] = [i + 'TTT' + str(Bio.Seq.Seq(i[B:]).reverse_complement())  for i in query['Sequence']]
            '''
            if B > 0:
                seqC = [str(Bio.Seq.Seq(i[:B]).reverse_complement()) for i in query['Sequence']]
            elif B < 0:
                seqC = [str(Bio.Seq.Seq(i[B:]).reverse_complement())  for i in query['Sequence']]
            if B!=0:
                queryC = query.copy()
                queryC['Sequence'] = seqC
                queryC['Strand'] = 'seq2' 
                query = pd.concat([query, queryC])
                query = query.sort_values(by=['name','pos1','pos2','Strand']) 
        else:
            query = df
            query['name'] = [i.split('/')[-1].split('.csv')[0] for i in df['filename']]

        # set concentrations
        query['Concentration (M)'] = 10e-9
        ref['Concentration (M)'] = 10e-9
        out = self.two_strand_analysis(query, ref)
        # fix names
        out['query'] = [i.split('_w')[0] for i in out['query']] 
        out['reference'] = [i.split('_w')[0] for i in out['reference']] 
        return out

    def two_strand_analysis(self, df1, df2, max_size=2):
        '''
        Run tube analysis on 2 strands
        '''
        logging.info('analyzing '+str(len(df1)) + 'x' + str(len(df2)))
        spec = nupack.SetSpec(max_size=max_size)
        df1['np_strand'] = [nupack.Strand(seq, name=name) for name, seq in df1[['Strand','Sequence']].values]
        df2['np_strand'] = [nupack.Strand(seq, name='seq0') for seq in df2['Sequence']]
        out = []
        for n1 in df1['name'].drop_duplicates().values:
            c = df1['name']==n1
            query_strand = [nps for nps in df1[c]['np_strand']]
            query_conc = [cnc for cnc in df1[c]['Concentration (M)']]
            [q_start, q_end] = df1[c][['pos1','pos2']].values[0]   
            for n2, ref_start, ref_end, ref_strand, ref_conc in df2[['name','pos1','pos2','np_strand','Concentration (M)']].values:
                logging.info('analyzing '+n1+' '+n2)
                # generate strands
                strands = [ref_strand]
                strands+= query_strand
                conc = [ref_conc]
                conc+= query_conc 
                # compute complexes
                cp_set = nupack.ComplexSet(strands=strands, complexes=spec)
                cp_res = nupack.complex_analysis(complexes=cp_set, model=self.model, compute=['pfunc'])
                # compute concentrations
                ts = {strands[i]:conc[i] for i in range(len(strands))}
                cc_res = nupack.complex_concentrations(tube=cp_set, data=cp_res, concentrations=ts)
                x = Analysis.get_dataframe(cp_res, cc_res)
                x['query'] = n1
                x['q_start'] = q_start 
                x['q_end'] = q_end
                x['reference'] = n2
                x['ref_start'] = ref_start
                x['ref_end'] = ref_end
                out.append(x)
        out = pd.concat(out)
        return out

    def cross_strand_analysis(self, df1, df2, max_size=2):
        '''
        Run tube analysis on 2 strands
        '''
        logging.info('analyzing '+str(len(df1)) + 'x' + str(len(df2)))
       
        # generate strands
        df1['np_strand'] = [nupack.Strand(seq, name=name) for name, seq in df1[['name','Sequence']].values]
        df2['np_strand'] = [nupack.Strand(seq, name=name) for name, seq in df2[['name','Sequence']].values]
       
        # generate complexes
        spec = Analysis.generate_complexes(df1, df2, csize) 
        cp_set = nupack.ComplexSet(strands=strands, complexes=spec)
        # run complex analysis
        cp_res = nupack.complex_analysis(complexes=cp_set, model=self.model, compute=['pfunc'])
        
        # run concentration analysis
        logging.info('computing concentrations '+str(len(df1))+'x'+str(len(df2)))
        out = []
        for n1, q_start, q_end in df1[['name','pos1','pos2']].values:
            for n2, ref_start, ref_end in df2[['name','pos1','pos2']].values:
                cc_res = nupack.complex_concentrations(tube=cp_set, data=cp_res, concentrations=ts)
                x = Analysis.get_dataframe(cc_res)
                x['query'] = n1
                x['q_start'] = q_start
                x['q_end'] = q_end
                x['reference'] = n2
                x['ref_start'] = ref_start
                x['ref_end'] = ref_end
                out.append(x)
        out = pd.concat(out)
        return out

    def off_target_score(self, df):
        '''
        Score off target
        '''
        df = df.groupby(['query','reference','q_start','q_end','Complex']).agg('sum').reset_index()
        c = df['query'] == df['reference']
        # compute on target and off target
        ontar = df[c].groupby(['query','q_start','q_end','Complex']).agg('sum').reset_index()
        offtar = df[~c].groupby(['query','q_start','q_end','Complex']).agg('sum').reset_index()
        col = ['query','q_start','q_end','Complex']
        ontar['reference'] = 'on_target'
        offtar['reference'] = 'off_target'
        offtar2 = offtar.copy()
        offtar2['Concentration (M)'] = 1/offtar2['Concentration (M)']
        offtar2['reference'] = 'inverse_off_target'
        df = pd.concat([df, ontar, offtar, offtar2])
        df = df.drop(columns=['ref_start','ref_end'])
        return df
 
    def ts25_spec(df):
        '''
        Filter designs satisfying the following reaction
        g1 + g2 -> g1g2
        g1g2 + t1 -> g1t1 + g2
        g1g2 + t1 -> g1 + t1g2
        '''
        c1 = df['Complex']=='g1[0]+g2[0]'
        c2 = df['Complex']=='g1[0]+wt1[0]'
        c2 = c2 | (df['Complex']=='g1[0]+5p')
        c2 = c2 | (df['Complex']=='g1[0]+3p')
        c3 = df['Complex']=='wt1[0]+g2[0]'
        c3 = c3 | (df['Complex']=='5p+g2[0]')
        c3 = c3 | (df['Complex']=='3p+g2[0]')

        s1 = df['Concentration (M)'] > 1e-9*9
        s2 = df['Concentration (M)'] < 1e-9*1

        t0 = df['tube']=='tube0'
        t1 = df['tube']=='tube1'

        f1 = df[s1 & c1 & t0]['filename']
        f2 = df[s2 & c1 & t1]['filename']

        df = df[df['filename'].isin(f1)]
        df = df[df['filename'].isin(f2)]
        return df

    def ts23_spec(df):
        '''
        Filter designs satisfying the following reaction
        g1 + g2 -> g1 + g2
        g1 + g2 + t1 -> g1g2t1
        '''
        c1 = df['Complex']=='g1[0]'
        c2 = df['Complex']=='g2[0]'
        '''
        c3 = df['Complex']=='wt1[0]+g1[0]+g2[0]'
        c3 = c3 | (df['Complex']=='5p+g1[0]+g2[0]')
        c3 = c3 | (df['Complex']=='3p+g1[0]+g2[0]')
        '''

        s1 = df['Concentration (M)'] > 1e-9*9
        s2 = df['Concentration (M)'] < 1e-9*1

        t0 = df['tube']=='tube0'
        t1 = df['tube']=='tube1'

        f1 = df[s1 & c1 & t0]['filename']
        f2 = df[s1 & c2 & t0]['filename']
        f3 = df[s2 & c1 & t1]['filename']
        f4 = df[s2 & c2 & t1]['filename']

        df = df[df['filename'].isin(f1)]
        df = df[df['filename'].isin(f2)]
        df = df[df['filename'].isin(f3)]
        df = df[df['filename'].isin(f4)]
        return df

    def ts45_spec(df):
        '''
        Filter designs satisfying the following reaction
        g1 + t1 -> g1t1
        '''
        c1 = df['Complex']=='g1[0]'
        c2 = df['Complex']=='g1[0]+wt1[0]'
        c2 = c2 | (df['Complex']=='wt1[0]+g1[0]')
        c2 = c2 | (df['Complex']=='g1[0]+5p')
        c2 = c2 | (df['Complex']=='g1[0]+3p')

        s1 = df['Concentration (M)'] > 1e-9*9
        s2 = df['Concentration (M)'] < 1e-9*1

        t0 = df['tube']=='tube0'
        t1 = df['tube']=='tube1'

        f1 = df[s1 & c1 & t0]['filename']
        f2 = df[s2 & c1 & t1]['filename']
        f3 = df[s1 & c2 & t1]['filename']

        df = df[df['filename'].isin(f1)]
        df = df[df['filename'].isin(f2)]
        df = df[df['filename'].isin(f3)]
        return df

    def filter_ts45(self, df):
        '''
        Sample function to filter for the N best designs for ts45
        '''
        # get trigger RNA
        trig = Oligos.get_tagged(df, ['t1*0*'])
        
        # set trigger window of 50bp on each side
        self.params['dim'] = [50]
        trig = self.get_ref_window(trig)
        trig['Strand'] = 'w'+trig['Strand']
        
        # merge dataframes
        col = ['Strand','Sequence','filename']
        dfin = pd.concat([trig[col], df[col]])
        dfin = dfin.sort_values(by=['filename','Strand'])
        
        # perform tube analysis for OFF state
        x = Oligos.get_tagged(dfin, ['*g1*0*'])
        self.params['dim'] = [10]
        self.params['maxsize'] = 2
        tube0 = self.on_target_analysis(x)
        # perform tube analysis for ON state
        x = Oligos.get_tagged(dfin, ['*g1*0*','wt1*0*'])
        self.params['dim'] = [10, 10]
        self.params['maxsize'] = 2
        tube1 = self.on_target_analysis(x)

        # filter for designs that meet test tube conditions
        tube0['tube'] = 'tube0'
        tube1['tube'] = 'tube1'
        tubes = pd.concat([tube0, tube1])
        # remove spec tags
        tubes['Complex'] = [i.replace('spec_','') for i in tubes['Complex']]
        tubes['Complex'] = [i.replace('chlor_','') for i in tubes['Complex']]
        tubes = Analysis.ts45_spec(tubes)
        
        N = len(tubes['filename'].drop_duplicates())
        total = len(df['filename'].drop_duplicates())
        logging.info('filtered '+str(N)+'/'+str(total))
        
        df = df[df['filename'].isin(tubes['filename'].drop_duplicates())]
        return df
        
    def plot_scatter(data, col, color='blue', scale=1):
        '''
        Scatter plot
        '''
        x1 = data['pos1'].astype(int)
        if 'pos2' in data.columns:
            x2 = data['pos2'].astype(int)
        else:
            L = np.array([len(i) for i in data['Sequence']])
            x2 = x1 + L
        y = data[col]
        x = data['pos1']
        dx = data['pos2'] - data['pos1']
        dy = [0]*len(x)
        plt.quiver(x, y, dx, dy, color=color, width=scale)

    def plot_line(data, col, color='blue', scale=1):
        '''
        Plot lines
        '''
        # get positions
        x1 = data['pos1'].astype(int).values
        if 'pos2' in data.columns:
            x2 = data['pos2'].astype(int).values
        else:
            L = np.array([len(i) for i in data['Sequence']])
            x2 = x1 + L
        y1 = data[col].values
        x, y = Analysis.get_window_average(x1, x2, y1)
        plt.plot(x, y*scale, color=color)

    def get_window_average(x1, x2, y1):
        '''
        Get a single line which is an average over window values
        '''
        # sum and average the values
        x = np.arange(min(x1),max(x2))
        c = np.zeros(len(x))
        y = np.zeros(len(x))
        for i in range(len(x2)):
            idx = (x >= x1[i]) & (x < x2[i])
            c[idx]+= 1
            y[idx]+= y1[i]
        y = y/c
        return x,y

    def get_2d_average(df, col):
        '''
        Get a 2d average over window values
        '''
        # generate meshgrid 
        y1 = min(df['q_start'])
        y1 = 0
        y2 = max(df['q_end'])
        y = np.arange(y1,y2)
        x1 = min(df['ref_start'])
        x1 = 0
        x2 = max(df['ref_end'])
        x = np.arange(x1,x2)
        pos = np.meshgrid(x,y)
        z = np.zeros((len(y), len(x)))
        c = np.zeros((len(y), len(x)))

        for val, y1, y2, x1, x2 in df[[col,'q_start','q_end','ref_start','ref_end']].values:
            idx = (pos[0] >= x1) & (pos[0] < x2) & (pos[1] >= y1) & (pos[1] < y2) 
            c[idx]+= 1
            z[idx]+= val
        c[c==0] = 1
        z = z/c
        return z

    def get_pos(self, df):
        df['q_start'] = [i.split('_')[-2] for i in df['query']]
        df['q_end'] = [i.split('_')[-1] for i in df['query']]
        df['ref_start'] = [i.split('_')[-2] for i in df['reference']]
        df['ref_end'] = [i.split('_')[-1] for i in df['reference']]
        df['query'] = [i.split('_')[0] for i in df['query']] 
        df['reference'] = [i.split('_')[0] for i in df['reference']] 
        return df
