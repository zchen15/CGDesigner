#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# library for designing cgRNAs with NUPACK

# computing libraries
import nupack
import numpy as np
import pandas as pd
import Bio
import Bio.Seq

# file and process management
import logging
import sys

# custom libraries
from .formula import Formula
from .misc import update_parameters
from .misc import read_csv

class Designer:
    '''
    Object holding functions to run test tube design
    '''
    codons = {'I':['ATT', 'ATC', 'ATA'],
            'L':['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
            'V':['GTT', 'GTC', 'GTA', 'GTG'],
            'F':['TTT', 'TTC'],
            'M':['ATG'],
            'C':['TGT', 'TGC'],
            'A':['GCT', 'GCC', 'GCA', 'GCG'],
            'G':['GGT', 'GGC', 'GGA', 'GGG'],
            'P':['CCT', 'CCC', 'CCA', 'CCG'],
            'T':['ACT', 'ACC', 'ACA', 'ACG'],
            'S':['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'Y':['TAT', 'TAC'],
            'W':['TGG'],
            'Q':['CAA', 'CAG'],
            'N':['AAT', 'AAC'],
            'H':['CAT', 'CAC'],
            'E':['GAA', 'GAG'],
            'D':['GAT', 'GAC'],
            'K':['AAA', 'AAG'],
            'R':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            '*':['TAA', 'TAG', 'TGA']}
    
    params = {'cache':2.0,
              'parallelize':False,
              'material':'rna',
              'temperature':37,
              'gin':'',
              'gout':'',
              'excludes':['g'*10,'c'*10,'a'*10,'t'*4],
              'enzymes':{'Esp3I':'CGTCTC',
                         'BsaI':'GGTCTC',
                         'BspMI':'ACCTGC',
                         'PaqCI':'CACCTGC'},
              'mtrig':[0]}

    def __init__(self, args=None):
        # load params
        self.params = update_parameters(args, self.params)
    
    def get_window(self, input_sequence, length, stride):
        '''
        Generate sequence windows of length and stride
        '''
        # remove spaces
        scan_window=[]
        for name, seq in input_sequence[['Strand','Sequence']].values:
            w1 = np.arange(0, len(seq), stride)
            w2 = w1 + length
            w2[w2 > len(seq)] = len(seq)
            w = [['w'+str(w1[i])+'_'+str(w2[i]), seq[w1[i]:w2[i]]] for i in range(len(w1))]
            scan_window+= [[name+'_'+n, s] for n, s in w]
        col = ['Strand','Sequence']
        scan_window = pd.DataFrame(scan_window, columns=col)
        scan_window['len'] = [len(i) for i in scan_window['Sequence']]
        scan_window['RNA'] = [self.isRNA(i) for i in scan_window['Sequence']]
        scan_window['AA'] = [self.isAminoAcid(seq) & ~rna for seq, rna in scan_window[['Sequence','RNA']].values]
        # prune out sequences not matching length
        for i in self.params['mtrig']:
            trigRNA = self.trigRNA[i]
            trigLen = np.sum([d.nt() for d in trigRNA[0]])
            c1 = scan_window['RNA'] & (scan_window['len'] < trigLen)
            scan_window = scan_window[~c1]
        if len(scan_window) == 0:
            logging.error('scan_window='+str(length)+' trigLen='+str(trigLen)) 
        return scan_window

    def get_constraints(self):
        '''
        Applies constraints to the cgRNA design
        '''
        # define soft constraints
        self.soft = []
        # define hard constraints
        self.hard = []
        # apply prevent patterns
        enzymes = [i for i in self.params['enzymes'].values()]
        enzymes+= [str(Bio.Seq.Seq(i).reverse_complement()) for i in enzymes]
        self.hard+= [nupack.Pattern(enzymes, scope=self.scope[j]) for j in range(len(self.scope))]
        homopolymers = self.params['excludes']
        self.hard+= [nupack.Pattern(homopolymers, scope=self.scope[j]) for j in range(len(self.scope))]

        # apply diversity constraint
        if len(self.gRNA_in) == 0:
            logging.info('Applying diversity constraint')
            self.hard+= [nupack.Diversity(word=8, types=4, scope=self.scope[j]) for j in range(len(self.scope))]

        # apply window constraint
        elif len(self.gRNA_in) > 0:
            for t in self.params['mtrig']:
                trigRNA = self.trigRNA[t]
                for j in range(len(trigRNA)):
                    name, seq = self.gRNA_in[['Strand','Sequence']].values[j]
                    seq = seq.upper().replace(' ','')
                    N = np.sum([d.nt() for d in trigRNA[j]])
                    # check if sequence is amino acids
                    try:
                        if self.isRNA(seq):
                            logging.info('RNA constraint for trig'+str(j)+str(t))
                            self.hard+=[nupack.Window(trigRNA[j], sources=[seq])]
                        elif self.isAminoAcid(seq):
                            logging.info('Amino Acid constraint for trig'+str(j))
                            self.hard+=[nupack.Library(trigRNA[j], self.get_AAlib(seq))]
                    except:
                        logging.error('sequence constraint '+name+' is invalid')
                        logging.error('sequence = '+seq)
                        logging.error('trigger length = '+str(N))
                        logging.error('mRNA length = '+str(len(seq)))
                        logging.error(str(self.params))
                        sys.exit(1)

    def isRNA(self, seq):
        '''
        Returns true if values are from RNA or DNA sequence
        '''
        return set(seq.upper()).issubset(set('ATGCU'))
 
    def isAminoAcid(self, seq):
        '''
        Return true if values from amino acid sequence
        '''
        return set(seq.upper()).issubset(set(self.codons))

    def get_AAlib(self, AAseq):
        '''
        Define a domain composed of amino acids
        '''
        # define a library of codons for each amino acid
        catalog = []
        for j in AAseq.upper():
            catalog.append(self.codons[j])
        return catalog

    def run(self):
        '''
        Run the nupack design
        '''
        # enable parallelism and change block cache size
        nupack.config.parallelism = self.params['parallelize']
        nupack.config.cache = self.params['cache']
        # specify model params
        material = self.params['material']
        temperature = self.params['temperature']
        fstop = self.params['fstop']
        self.model = nupack.Model(material=material, celsius=temperature)
        self.options = nupack.DesignOptions(f_stop=fstop, wobble_mutations=True) # always a number in (0,1)
        
        # set guide RNA output targets
        df = read_csv(self.params['gout'])
        spec = Formula(N=self.params['Northo'], targets=df)

        # set test tube formula
        tube_spec = getattr(spec, self.params['spec'])
        tube_spec(self.params['dim'])
        formula = spec.get_formula() # extra data from formula
        self.tubes = formula['tubes']
        self.scope = formula['scopes']
        self.trigRNA = formula['trigger']

        # output filename
        outname = self.params['outname'].split('.')[0]
        outname+='_dim'+''.join([str(i)+'_' for i in self.params['dim']])
        outname = outname[:-1]

        # load input trigger sequence constraints
        df = read_csv(self.params['gin'])
        if 'scan' in self.params.keys():
            window, stride = self.params['scan']
            df = self.get_window(df, window, stride)
            
        designs = []
        # run single design
        if type(df)==list and df==[]:
            self.gRNA_in = []
            self.get_constraints()
            d = nupack.tube_design(self.tubes, soft_constraints=self.soft, hard_constraints=self.hard, options=self.options, model=self.model)
            fname = outname+'.spec'
            logging.info('exporting design to '+fname)
            d.save(fname)

        # run multiple designs
        for i in range(len(df)):
            # apply constraints
            seqname, sequence = df[['Strand','Sequence']].values[i]
            self.gRNA_in = pd.DataFrame([[seqname, sequence]], columns=['Strand','Sequence'])
            self.get_constraints()
            # generate new designs
            d = nupack.tube_design(self.tubes, soft_constraints=self.soft, hard_constraints=self.hard, options=self.options, model=self.model)
            fname = outname+'_'+seqname+'.spec'
            logging.info('exporting design to '+fname)
            d.save(fname)


