#!/usr/bin/python3 bash

import subprocess
import numpy as np

def subcall(cmd):
    for c in cmd:
        print(c)
        subprocess.call(c, shell=True)

# Generates a set of designs for the inhibited split terminator switch targeting a set of mRNA sequences 
dimA = np.arange(5,30,1) # dimensions for toehold A
dimB = np.arange(3,30,1) # dimensions for duplex B
cmd = []
for a in dimA:
    for b in dimB:
        cmd+= ['python3 main.py design -material rna -s inhibited_split_terminator_switch -gin sequences/mRNA.csv -o ts37_mRNA_'+str(len(cmd))+' -d 20 5 10 0 '+str(a)+' '+str(b)+' 10 6 -fstop 0.05 -scan 200 50 -pk']
#subcall(cmd)

# Generates a set of designs for the reverse toehold switch targeting a set of mRNA sequences
dimA = [15,30]    # dimension for toehold A
dimB = [10,15,20] # dimension for duplex B
cmd = []
for a in dimA:
    for b in dimB:
        cmd+= ['python3 main.py design -material rna -s reverse_toehold_switch -gin sequences/PAX7.csv -o ts45_mRNA_'+str(len(cmd))+' -d '+str(b)+' '+str(a)+' 6 3 6 3 1 6 -fstop 0.05 -scan 100 10 -pk']
subcall(cmd)

