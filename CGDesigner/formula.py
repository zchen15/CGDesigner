#!/usr/bin/env python
# Zhewei Chen (zchen@caltech.edu)
# Objects holding information about cgRNA test tube formulas

import nupack
import numpy as np

class Formula:
    '''
    This class object holds data on test formulations for cgRNAs
    '''
    params = {'parts':{'t1Tm':{'sequence':'cgcaaaaaaccccgcttcggcggggttttttcgc',
                               'DU':'U1D2(D12(U4)U1)'},
                       't2Tm':{'sequence':'ctgtaacagagcattagcgcaaggtgatttttgtcttcttgcgctaatttttt',
                               'DU':'U7D4(U1D12(U1D3(U4)U1)U1)'},
                       'g1Tm':{'sequence':'cgccgcaaaccccgcccctgacagggcggggtttcgccgc',
                               'DU':'U1D2(U1D2(D11(U6)U1)U1)'},
                       'g2Tm':{'sequence':'aaaaaaaaaccccgcccctgacagggcggggtttttttt',
                               'DU':'U1D16(U6)'}}}

    def __init__(self, N=1, targets=[]):
        '''
        Initializes the design formula
        N = number of orthogonal cgRNAs to design
        target = dna targeting sequence
        '''
        # set input targets and number of cgRNA to design
        if targets!=[]:
            self.N = len(targets)
            self.target = [nupack.Domain(seq, name=name) for name, seq in targets[['name','Sequence']].values]
        else:
            self.N=N
            self.target = [nupack.Domain('AACTTTCAGTTTAGCGGTCT', name='xRFP')]*self.N # set xRFP as default dna target

        # Cas9 domains
        self.cas9h = nupack.Domain('gttttAGAgctaGAAAtagcAAGTTaaaatAAggCTAGTccG', name='cas9h')
        self.DU_cas9h = 'U1D6(U1 D4(U4) U3)U1 U2D2(U5)U1'
        self.cas9_a = nupack.Domain('gttttAGA', name='cas9_a')
        self.cas9_b = nupack.Domain('AAGTTaaaat', name='cas9_b')
        self.cas9_c = nupack.Domain('AAggCTAGTccG', name='cas9_c')
        self.DU_cas9_c = 'U2D2(U5)U1'
        self.cas9ha = nupack.Domain('gttttAGAgcta', name='cas9ha')
        self.cas9hb = nupack.Domain('tagcAAGTTaaaatAAggCTAGTccG', name='cas9hb')
        self.DU_cas9ha = 'U1D6(U1 D4('
        self.DU_cas9hb = ') U3)U1 U2D2(U5)U1'
        # S.Pyogene terminator
        self.spyTerm = nupack.Domain('ttatcaACTTgaaaAAGTgGCACCGagtCGGTGCttttttt', name='spyTerm')
        self.DU_spyTerm = 'U6 D4(U4) U1D6(U3)U7' 
        # Terminators
        self.term2 = nupack.Domain('ggcaccgagtcggtgcttttttt', name='Term2')
        self.DU_term2 = 'U1D6(U3)U7'
        self.g1Tm = nupack.Domain('cgccgcaaaccccgcccctgacagggcggggtttcgccgc', name='g1Tm')
        self.DU_g1Tm = 'U1D2(U1D2(D11(U6)U1)U1)'
        self.g2Tm = nupack.Domain('aaaaaaaaaccccgcccctgacagggcggggtttttttt', name='g2Tm')
        self.DU_g2Tm = 'U1D16(U6)'
        self.t1Tm = nupack.Domain('cgcaaaaaaccccgcttcggcggggttttttcgc', name='t1Tm')
        self.DU_t1Tm = 'U1D2(D12(U4)U1)'
        self.t2Tm = nupack.Domain('ctgtaacagagcattagcgcaaggtgatttttgtcttcttgcgctaatttttt', name='t2Tm')
        self.DU_t2Tm = 'U7D4(U1D12(U1D3(U4)U1)U1)'

    def get_formula(self):
        '''
        Get design formula
        return tubes and constraints related to each design
        '''
        data = {'tubes':self.tubes,
                'trigger':self.trigRNA,
                'scopes':self.scope,
                'weights':self.weights}
        return data
        
    def rna_trig(self, dimensions=[10]):
        '''
        Sets the formula to design an orthogonal set of RNA triggers
        '''
        # variable domains
        [A] = dimensions
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        
        # define domain scopes for soft and hard constraints
        self.scope = [[dA[j]] for j in range(self.N)]

        # define domain scopes for trigger RNA
        self.trigRNA = [[[dA[j]] for j in range(self.N)]]

        # strands
        s1 = [nupack.TargetStrand([dA[i]], name=['sF_'+str(i)]) for i in range(self.N)]
        s2 = [nupack.TargetStrand([~dA[i]], name=['sR_'+str(i)]) for i in range(self.N)]

        # target complexes
        c3 = [nupack.TargetComplex([s1[i]], structure='U'+str(dA[i].nt()), name=['c1',i]) for i in range(self.N)]
        c4 = [nupack.TargetComplex([s2[i]], structure='U'+str(dA[i].nt()), name=['c2',i]) for i in range(self.N)]
        
        self.tubes = [] 
        # crosstalk
        for i in range(self.N):
            for j in range(self.N):
                if i!=j:
                    cross = [nupack.TargetTube({c3[i]:500e-9, c4[j]:500e-9}, nupack.SetSpec(max_size=2), name=['cross',i,j])]
                    self.tubes+= cross
        # weights
        weights = nupack.Weights(self.tubes)
        for i in range(self.N):
            weights[dA[i], :, :, :] = 1
        self.weights = weights

    def catalytic_two_hairpin_switch(self, dimensions=[10,10,10,10,10,6]):
        '''
        Sets design specifications for ts12 (catalytic two hairpin switch)
        OFF -> ON trigger logic
        g1, g2 are guide RNAs split at the terminator stem
        g1 + g2 -> g1 + g2
        g1 + t1 -> t1g1 + g2 -> t1g1g2 + g1 -> t1g1g2g1 + g2 -> t1g1g2g1g2 ...
 
        t1g1g2 is an active guide RNA
        g1,g2 are sequestered by a hairpin
        '''
        # unpack dimensions
        [A,B,C,L] = dimensions
        # fixed domains
        target = self.target
        DU_cas9h = self.DU_cas9h
        cas9h = self.cas9h 
        term2 = self.term2
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        DU_term2 = self.DU_term2
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        # specify strands
        Sg1 = [nupack.TargetStrand([~dA[j], ~dB[j], target[j], cas9h, dL[j], dC[j], dB[j], g1Tm], name=['spec_g1',j]) for j in range(self.N)]
        Sg2 = [nupack.TargetStrand([dB[j], dA[j], ~dB[j], ~dC[j], term2, g2Tm], name=['chlor_g2',j]) for j in range(self.N)]
        St1 = [nupack.TargetStrand([dB[j], dA[j], t1Tm], name=['t1',j]) for j in range(self.N)]
        St2 = [nupack.TargetStrand([dC[j], dB[j], t1Tm], name=['t2',j]) for j in range(self.N)]
 
        Cg1 = [nupack.TargetComplex([Sg1[j]],
                structure='U'+str(dA[j].nt())+'D'+str(dB[j].nt())+'(U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt()+dC[j].nt())+')'+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Cg2 = [nupack.TargetComplex([Sg2[j]],
                structure='D'+str(dB[j].nt())+'(U'+str(dA[j].nt())+')U'+str(dC[j].nt())+DU_term2+DU_g2Tm,
                name=['Cg2', j]) for j in range(self.N)]
        Ct1 = [nupack.TargetComplex([St1[j]],
                structure='U'+str(dB[j].nt()+dA[j].nt())+DU_t1Tm,
                name=['Ct1',j]) for j in range(self.N)]

        Ct2 = [nupack.TargetComplex([St2[j]],
                structure='U'+str(dC[j].nt()+dB[j].nt())+DU_t1Tm,
                name=['Ct2',j]) for j in range(self.N)]

        Ct1g1 = [nupack.TargetComplex([St1[j], Sg1[j]],
                structure='D'+str(dB[j].nt()+dA[j].nt())+'('+DU_t1Tm+'+) U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt()+dC[j].nt()+dB[j].nt())+DU_g1Tm,
                name=['Ct1g1', j]) for j in range(self.N)]

        Ct2g2 = [nupack.TargetComplex([St2[j], Sg2[j]],
                structure='D'+str(dC[j].nt()+dB[j].nt())+'('+DU_t1Tm+'+ U'+str(dB[j].nt()+dA[j].nt())+')'+DU_term2+DU_g2Tm,
                name=['Ct2g2', j]) for j in range(self.N)]

        Ct1g1g2 = [nupack.TargetComplex([St1[j], Sg1[j], Sg2[j]],
                structure='D'+str(dB[j].nt()+dA[j].nt())+'('+DU_t1Tm+'+) U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+'D'+str(dC[j].nt()+dB[j].nt())+'('+DU_g1Tm+'+U'+str(dB[j].nt()+dA[j].nt())+')'+DU_term2+DU_g2Tm,
                name=['Ct1g1g2', j]) for j in range(self.N)]
        
        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + g2 -> g1 + g2
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=2), name=['step0',j]))
            # g1 + t1 -> t1g1
            steps.append(nupack.TargetTube({Ct1[j]:1e-8, Cg1[j]:1e-8}, nupack.SetSpec(max_size=2, exclude=[Ct1g1[j]]), name=['step1a',j]))
            # g2 + t2 -> t2g2
            steps.append(nupack.TargetTube({Ct2[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=2, exclude=[Ct2g2[j]]), name=['step1b',j]))
            # t1g1 + g2 -> t1g1g2
            steps.append(nupack.TargetTube({Ct1g1[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=3, exclude=[Ct1g1g2[j]]), name=['step2',j]))
            # t1g1g2 -> t1g1g2
            steps.append(nupack.TargetTube({Ct1g1g2[j]:1e-8}, nupack.SetSpec(max_size=3), name=['step3',j]))
        self.tubes = steps
        
        # set the weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights
        
        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dC[j], dB[j]] for j in range(self.N)]
        self.scope+= [[dB[j], dA[j], ~dB[j], ~dC[j]] for j in range(self.N)]
        
        # define trigRNA domains
        self.trigRNA = [[[dB[j], dA[j]] for j in range(self.N)],
                   [[dC[j], dB[j]] for j in range(self.N)]]

    def single_hairpin_switch(self, dimensions=[5,20,10,6]):
        '''
        Sets design specifications for ts18b (single hairpin switch)
        OFF -> ON trigger logic
        g1 and g2 are split at the terminator
        g1 + t1 -> t1g1 + g2 -> + t1g1g2
        '''
        # unpack dimensions
        [A,B,C,E,F,G,L] = dimensions
        
        # fixed domains
        target = self.target
        cas9h = self.cas9h
        DU_cas9h = self.DU_cas9h
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        t2Tm = self.t2Tm
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        DU_t2Tm = self.DU_t2Tm
        
        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        
        # specify strands
        term2 = [[dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]

        Sg1 = [nupack.TargetStrand([~dC[j], ~dB[j], target[j], cas9h, dL[j], dA[j], dB[j], g1Tm], name=['spec_g1',j]) for j in range(self.N)]
        Sg2 = [nupack.TargetStrand([~dA[j]] + term2[j] + [g2Tm], name=['chlor_g2',j]) for j in range(self.N)]
        St1 = [nupack.TargetStrand([dB[j], dC[j], t1Tm], name=['t1',j]) for j in range(self.N)]

        Cg1 = [nupack.TargetComplex([Sg1[j]],
                structure='U'+str(dC[j].nt())+'D'+str(dB[j].nt())+'(U'+str(target[j].nt())+DU_cas9h+' U'+str(dL[j].nt()+dA[j].nt())+')'+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Cg2 = [nupack.TargetComplex([Sg2[j]],
                structure='U'+str(dA[j].nt())+DU_term2[j]+DU_g2Tm,
                name=['Cg2', j]) for j in range(self.N)]
        Ct1 = [nupack.TargetComplex([St1[j]],
                structure='U'+str(dB[j].nt()+dC[j].nt())+DU_t1Tm,
                name=['Ct1',j]) for j in range(self.N)]
        Ct1g1 = [nupack.TargetComplex([St1[j], Sg1[j]],
                structure='D'+str(dC[j].nt()+dB[j].nt())+'('+DU_t1Tm+'+) U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt()+dA[j].nt()+dB[j].nt())+DU_g1Tm,
                name=['Ct1g1',j]) for j in range(self.N)]
        Ct1g1g2 = [nupack.TargetComplex([St1[j], Sg1[j], Sg2[j]],
                structure='D'+str(dC[j].nt()+dB[j].nt())+'('+DU_t1Tm+'+) U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+'D'+str(dA[j].nt())+'(U'+str(dB[j].nt())+DU_g1Tm+'+)'+DU_term2[j]+DU_g2Tm,
                name=['Ct1g1g2',j]) for j in range(self.N)]
 
        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + g2 -> g1 + g2
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=2), name=['step0',j]))
            # t1 + g1 -> g1t1
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Ct1[j]:1e-8}, nupack.SetSpec(max_size=2, exclude=[Ct1g1[j]]), name=['step1',j]))
            # t1g1 + g2 -> t1g1g2
            steps.append(nupack.TargetTube({Ct1g1[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=3, exclude=[Ct1g1g2[j]]), name=['step2',j]))
            # t1g1g2
            steps.append(nupack.TargetTube({Ct1g1g2[j]:1e-8}, nupack.SetSpec(max_size=3), name=['step3',j]))
        self.tubes = steps
        
        # set the weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights

        # define domain scopes for soft and hard constraints
        self.scope = [[dB[j], dC[j]] for j in range(self.N)]
        self.scope+= [[dL[j], dA[j], dB[j]] for j in range(self.N)]
        self.scope+= [[~dA[j]] + term2[j] for j in range(self.N)]

        # define trigRNA domains
        self.trigRNA = [[[dB[j], dC[j]] for j in range(self.N)]]

    def two_hairpin_switch(self, dimensions=[10,10,5,5,5,5,6,1,7,3,6]):
        '''
        Sets design specifications for ts23b (two hairpin switch)
        OFF -> ON trigger logic
        g1, g2 are guide RNAs split at the terminator stem
        g1 + g2 -> g1 + g2
        g1 + t1 -> t1g1 + g2 -> t1g1g2
 
        t1g1g2 is an active guide RNA
        g1,g2 are sequestered by a hairpin
        '''
        # unpack dimensions
        [A,B,C,x,y,z,L,D,E,F,G] = dimensions

        # fixed domains
        target = self.target
        DU_cas9h = self.DU_cas9h
        cas9h = self.cas9h 
        term2 = self.term2
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        
        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dX = [nupack.Domain('N'*x, name=['X',j]) for j in range(self.N)]
        dY = [nupack.Domain('N'*y, name=['Y',j]) for j in range(self.N)]
        dZ = [nupack.Domain('N'*z, name=['Z',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        dD = [nupack.Domain('N'*D, name=['D',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        term2 = [[dD[j], dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['U'+str(dD[j].nt())+'D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]

        # specify strands
        Sg1 = [nupack.TargetStrand([~dA[j], ~dX[j], ~dB[j], target[j], cas9h, dL[j], dC[j], dY[j], dB[j], dX[j], g1Tm], name=['spec_g1',j]) for j in range(self.N)]
        Sg2 = [nupack.TargetStrand([dY[j], dB[j], dZ[j], ~dB[j], ~dY[j], ~dC[j]]+term2[j]+[g2Tm], name=['chlor_g2',j]) for j in range(self.N)]
        St1 = [nupack.TargetStrand([dB[j], dX[j], dA[j], t1Tm], name=['t1',j]) for j in range(self.N)]
 
        Cg1 = [nupack.TargetComplex([Sg1[j]],
                structure='U'+str(dA[j].nt())+'D'+str(dX[j].nt()+dB[j].nt())+'(U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt()+dC[j].nt()+dY[j].nt())+')'+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Cg2 = [nupack.TargetComplex([Sg2[j]],
                structure='D'+str(dY[j].nt()+dB[j].nt())+'(U'+str(dZ[j].nt())+')U'+str(dC[j].nt())+DU_term2[j]+DU_g2Tm,
                name=['Cg2', j]) for j in range(self.N)]
        Ct1 = [nupack.TargetComplex([St1[j]],
                structure='U'+str(dB[j].nt()+dX[j].nt()+dA[j].nt())+DU_t1Tm,
                name=['Ct1',j]) for j in range(self.N)]

        Ct1g1 = [nupack.TargetComplex([St1[j], Sg1[j]],
                structure='D'+str(dB[j].nt()+dX[j].nt()+dA[j].nt())+'('+DU_t1Tm+'+) U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt()+dC[j].nt()+dY[j].nt()+dB[j].nt()+dX[j].nt())+DU_g1Tm,
                name=['Ct1g1', j]) for j in range(self.N)]
        Ct1g1g2 = [nupack.TargetComplex([St1[j], Sg1[j], Sg2[j]],
                structure='D'+str(dB[j].nt()+dX[j].nt()+dA[j].nt())+'('+DU_t1Tm+'+) U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+'D'+str(dC[j].nt()+dY[j].nt()+dB[j].nt())+'(U'+str(dX[j].nt())+DU_g1Tm+'+U'+str(dY[j].nt()+dB[j].nt()+dZ[j].nt())+')'+DU_term2[j]+DU_g2Tm,
                name=['Ct1g1g2', j]) for j in range(self.N)]
        
        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + g2 -> g1 + g2
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=2), name=['step0',j]))
            # g1 + t1 -> t1g1
            steps.append(nupack.TargetTube({Ct1[j]:1e-8, Cg1[j]:1e-8}, nupack.SetSpec(max_size=2, exclude=[Ct1g1[j]]), name=['step1',j]))
            # t1g1 + g2 -> t1g1g2
            steps.append(nupack.TargetTube({Ct1g1[j]:1e-8, Cg2[j]:1e-8}, nupack.SetSpec(max_size=3, exclude=[Ct1g1g2[j]]), name=['step2',j]))
            # t1g1g2 -> t1g1g2
            steps.append(nupack.TargetTube({Ct1g1g2[j]:1e-8}, nupack.SetSpec(max_size=3), name=['step3',j]))
        self.tubes = steps
        
        # set the weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights
        
        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dC[j], dY[j], dB[j], dX[j]] for j in range(self.N)]
        self.scope+= [[dY[j], dB[j], dZ[j], ~dB[j], ~dY[j], ~dC[j]]+term2[j] for j in range(self.N)]

        # define trigRNA domains
        self.trigRNA = [[[dB[j], dX[j], dA[j]] for j in range(self.N)]]

    def reverse_split_terminator_switch(self, dimensions=[20,10,10,1,7,3,7,6]):
        '''
        Sets design specifications for ts25 (reverse split terminator switch)
        ON -> OFF trigger logic
        g1 + g2 -> g1g2 + t1 + t2 -> g1t1 + t2g2
        g1g2 is an active guide RNA
        t1 or t2 inactivates the guide RNA
        guide RNA is split at terminator
        '''
        # unpack dimensions
        [A,B,C,D,E,F,G,L] = dimensions

        # fixed domains
        target = self.target
        cas9h = self.cas9h
        DU_cas9h = self.DU_cas9h
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        t2Tm = self.t2Tm
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        DU_t2Tm = self.DU_t2Tm

        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dD = [nupack.Domain('N'*D, name=['D',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        term2 = [[dD[j], dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['U'+str(dD[j].nt())+'D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]
 
        # specify strands
        Sg1 = [nupack.TargetStrand([target[j], cas9h, dL[j], dA[j], dB[j], g1Tm], name=['chlor_g1',j]) for j in range(self.N)]
        Sg2 = [nupack.TargetStrand([~dC[j], ~dA[j]]+term2[j]+[g2Tm], name=['spec_g2',j]) for j in range(self.N)]
        St1 = [nupack.TargetStrand([~dB[j], ~dA[j], ~dL[j], t1Tm], name=['t1',j]) for j in range(self.N)]
        St2 = [nupack.TargetStrand([dA[j], dC[j], t2Tm], name=['t2',j]) for j in range(self.N)]
        
        # target complexes
        Cg1 = [nupack.TargetComplex([Sg1[j]],
                structure='U'+str(target[j].nt())+DU_cas9h+' U'+str(dL[j].nt())+' U'+str(dA[j].nt())+' U'+str(dB[j].nt())+' '+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Cg2 = [nupack.TargetComplex([Sg2[j]],
                structure='U'+str(dC[j].nt())+' U'+str(dA[j].nt())+DU_term2[j]+DU_g2Tm,
                name=['Cg2', j]) for j in range(self.N)]
        Cg1g2 = [nupack.TargetComplex([Sg1[j], Sg2[j]],
                structure='U'+str(target[j].nt())+DU_cas9h+' U'+str(dL[j].nt())+' D'+str(dA[j].nt())+'(U'+str(dB[j].nt())+' '+DU_g1Tm+'+ U'+str(dC[j].nt())+') '+DU_term2[j]+DU_g2Tm,
                name=['Cg1g2',j]) for j in range(self.N)]
        Ct1 = [nupack.TargetComplex([St1[j]],
                structure='U'+str(dB[j].nt())+'U'+str(dA[j].nt())+'U'+str(dL[j].nt())+DU_t1Tm,
                name=['Ct1',j]) for j in range(self.N)]
        Ct2 = [nupack.TargetComplex([St2[j]],
                structure='U'+str(dA[j].nt())+'U'+str(dC[j].nt())+DU_t2Tm,
                name=['Ct2',j]) for j in range(self.N)]
        Cg1t1 = [nupack.TargetComplex([Sg1[j], St1[j]],
                structure='U'+str(target[j].nt())+DU_cas9h+' D'+str(dL[j].nt())+'D'+str(dA[j].nt())+'D'+str(dB[j].nt())+'('+DU_g1Tm+'+) '+DU_t1Tm,
                name=['Cg1t1',j]) for j in range(self.N)]
        Ct2g2 = [nupack.TargetComplex([St2[j], Sg2[j]],
                structure='D'+str(dA[j].nt())+'D'+str(dC[j].nt())+'('+DU_t2Tm+'+) '+DU_term2[j]+DU_g2Tm,
                name=['Ct2g2',j]) for j in range(self.N)]

        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + g2 -> g1g2
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Cg2[j]:1e-8},
                          nupack.SetSpec(max_size=2, exclude=[Cg1g2[j]]),
                          name=['step0',j]))
            # g1g2 + t1 -> g1t1 + g2
            steps.append(nupack.TargetTube({Cg1g2[j]:1e-8, Ct1[j]:1e-8},
                          nupack.SetSpec(max_size=2, exclude=[Cg1t1[j], Cg2[j]]),
                          name=['step1a',j]))
            # g1t1 + g2 -> g1t1 + g2
            steps.append(nupack.TargetTube({Cg1t1[j]:1e-8, Cg2[j]:1e-8},
                          nupack.SetSpec(max_size=2),
                          name=['step2a',j]))
            # g1g2 + t2 -> g1 + t2g2
            steps.append(nupack.TargetTube({Cg1g2[j]:1e-8, Ct2[j]:1e-8},
                          nupack.SetSpec(max_size=2, exclude=[Cg1[j], Ct2g2[j]]),
                          name=['step1b',j]))
            # g1 + t2g2 -> g1 + t2g2
            steps.append(nupack.TargetTube({Ct2g2[j]:1e-8, Cg1[j]:1e-8},
                          nupack.SetSpec(max_size=2),
                          name=['step2b',j]))
        self.tubes = steps

        # design weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights

        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dA[j], dB[j]] for j in range(self.N)]
        self.scope+= [[~dC[j], ~dA[j], dD[j], dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        
        # define domain scopes for trigger RNA
        self.trigRNA = [[[~dB[j], ~dA[j], ~dL[j]] for j in range(self.N)],
                        [[dA[j], dC[j]] for j in range(self.N)]]

    def split_guideRNA(self, dimensions=[20,10,10,1,7,3,7,6]):
        '''
        Sets design specifications for ts26 (split terminator switch)
        ON -> OFF trigger logic
        g1 + g2 -> g1g2 + t1 + t2 -> g1t1 + t2g2
        g1g2 is an active guide RNA
        guide RNA is split at terminator
        '''
        # unpack dimensions
        [A,B,C,D,E,F,G,L] = dimensions
        
        # fixed domains
        target = self.target
        cas9h = self.cas9h
        DU_cas9h = self.DU_cas9h
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        t2Tm = self.t2Tm
        DU_term2 = self.DU_term2
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        DU_t2Tm = self.DU_t2Tm

        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dD = [nupack.Domain('N'*D, name=['D',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        term2 = [[dD[j], dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['U'+str(dD[j].nt())+'D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]
 
        # specify strands
        Sg1 = [nupack.TargetStrand([target[j], cas9h, dL[j], dA[j], dB[j], g1Tm], name=['chlor_g1',j]) for j in range(self.N)]
        Sg2 = [nupack.TargetStrand([~dC[j], ~dA[j]] + term2[j] + [g2Tm], name=['t1',j]) for j in range(self.N)]
        
        # target complexes
        Cg1 = [nupack.TargetComplex([Sg1[j]],
                structure='U'+str(target[j].nt())+DU_cas9h+' U'+str(dL[j].nt())+' U'+str(dA[j].nt())+' U'+str(dB[j].nt())+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Cg2 = [nupack.TargetComplex([Sg2[j]],
                structure='U'+str(dC[j].nt())+' U'+str(dA[j].nt())+DU_term2[j]+DU_g2Tm,
                name=['Cg2', j]) for j in range(self.N)]
        Cg1g2 = [nupack.TargetComplex([Sg1[j], Sg2[j]],
                structure='U'+str(target[j].nt())+DU_cas9h+' U'+str(dL[j].nt())+' D'+str(dA[j].nt())+'(U'+str(dB[j].nt())+DU_g1Tm+'+ U'+str(dC[j].nt())+') '+DU_term2[j]+DU_g2Tm,
                name=['Cg1g2',j]) for j in range(self.N)]

        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + g2 -> g1g2
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Cg2[j]:1e-8},
                          nupack.SetSpec(max_size=2, exclude=[Cg1g2[j]]),
                          name=['step0',j]))
            # g1g2 -> g1g2
            steps.append(nupack.TargetTube({Cg1g2[j]:1e-8},
                          nupack.SetSpec(max_size=2),
                          name=['step1',j]))
        self.tubes = steps
        
        # design weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights

        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dA[j], dB[j]] for j in range(self.N)]
        
        # define domain scopes for trigger RNA
        self.trigRNA = [[[~dC[j], ~dA[j]] + term2[j] for j in range(self.N)]]

    def terminator_switch(self, dimensions=[10,20,6,0,7,3,1]):
        '''
        Sets design specifications for ts32 (terminator switch)
        ON -> OFF trigger logic
        g1 is a terminator switch guide RNA
        g1 + t1 -> t1g1
        g1 is sequestered by t1
        '''
        # unpack dimensions
        [A,B,L,D,E,F,G] = dimensions

        # fixed domains
        target = self.target
        DU_cas9h = self.DU_cas9h
        cas9h = self.cas9h 
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        
        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        dD = [nupack.Domain('N'*D, name=['D',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        term2 = [[dD[j], dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['U'+str(dD[j].nt())+'D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]
        
        # specify strands
        Sg1 = [nupack.TargetStrand([target[j], cas9h, dL[j], dA[j], dB[j], ~dA[j]] + term2[j] + [g1Tm], name=['spec_g1',j]) for j in range(self.N)]
        St1 = [nupack.TargetStrand([~dB[j], ~dA[j], ~dL[j], t1Tm], name=['t1',j]) for j in range(self.N)]

        Cg1 = [nupack.TargetComplex([Sg1[j]],
            structure='U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+'D'+str(dA[j].nt())+'(U'+str(dB[j].nt())+')'+DU_term2[j]+DU_g1Tm,
            name=['Cg1', j]) for j in range(self.N)]
        Ct1 = [nupack.TargetComplex([St1[j]],
            structure='U'+str(dB[j].nt()+dA[j].nt()+dL[j].nt())+DU_t1Tm,
            name=['Ct1',j]) for j in range(self.N)]
        Cg1t1 = [nupack.TargetComplex([Sg1[j], St1[j]], structure='U'+str(target[j].nt())+DU_cas9h+'D'+str(dL[j].nt()+dA[j].nt()+dB[j].nt())+'(U'+str(dA[j].nt())+DU_term2[j]+DU_g1Tm+'+)'+DU_t1Tm,
            name=['Ct1g1', j]) for j in range(self.N)]
        
        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + t1 -> g1t1
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Ct1[j]:1e-8}, nupack.SetSpec(max_size=2, exclude=[Cg1t1[j]]), name=['step0',j]))
            # g1t1 -> g1t1
            steps.append(nupack.TargetTube({Cg1t1[j]:1e-8}, nupack.SetSpec(max_size=2), name=['step1',j]))
        self.tubes = steps

        # set the weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights
        
        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dA[j], dB[j]] for j in range(self.N)]
        self.scope+= [term2[j] for j in range(self.N)]

        # define trigRNA domains
        self.trigRNA = [[[dA[j], ~dB[j]] for j in range(self.N)]]

    def inhibited_split_terminator_switch(self, dimensions=[20,10,10,1,7,3,7,6]):
        '''
        Sets design specifications for ts38 (inhibited split terminator switch)
        OFF -> ON trigger logic
        g1 + g2 -> g1g2
        g1g2 is an active guide RNA
        guide RNA is split at terminator
        '''
        # unpack dimensions
        [A,B,C,D,E,F,G,L] = dimensions

        # fixed domains
        target = self.target
        cas9h = self.cas9h
        DU_cas9h = self.DU_cas9h
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        t2Tm = self.t2Tm
        DU_term2 = self.DU_term2
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        DU_t2Tm = self.DU_t2Tm

        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dD = [nupack.Domain('N'*D, name=['D',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        term2 = [[dD[j], dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['U'+str(dD[j].nt())+'D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]
 
        # specify strands
        Sg1 = [nupack.TargetStrand([~dA[j], target[j], cas9h, dL[j], dA[j], dB[j], g1Tm], name=['chlor_g1',j]) for j in range(self.N)]
        Sg2 = [nupack.TargetStrand([~dC[j], ~dB[j], ~dA[j]] + term2[j] + [g2Tm], name=['t1',j]) for j in range(self.N)]
        
        # target complexes
        Cg1 = [nupack.TargetComplex([Sg1[j]],
                structure='D'+str(dA[j].nt())+'(U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+')U'+str(dB[j].nt())+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Cg2 = [nupack.TargetComplex([Sg2[j]],
                structure='U'+str(dC[j].nt()+dB[j].nt()+dA[j].nt())+DU_term2[j]+DU_g2Tm,
                name=['Cg2', j]) for j in range(self.N)]
        Cg1g2 = [nupack.TargetComplex([Sg1[j], Sg2[j]],
                structure='U'+str(dA[j].nt()+target[j].nt())+DU_cas9h+' U'+str(dL[j].nt())+'D'+str(dA[j].nt()+dB[j].nt())+'('+DU_g1Tm+'+ U'+str(dC[j].nt())+')'+DU_term2[j]+DU_g2Tm,
                name=['Cg1g2',j]) for j in range(self.N)]

        # test tubes
        steps = []
        for j in range(self.N):
            # g1 + g2 -> g1g2
            steps.append(nupack.TargetTube({Cg1[j]:1e-8, Cg2[j]:1e-8},
                          nupack.SetSpec(max_size=2, exclude=[Cg1g2[j]]),
                          name=['step0',j]))
            # g1g2 -> g1g2
            steps.append(nupack.TargetTube({Cg1g2[j]:1e-8},
                          nupack.SetSpec(max_size=2),
                          name=['step1',j]))
        self.tubes = steps
        
        # design weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights

        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dA[j], dB[j]] for j in range(self.N)]
        
        # define domain scopes for trigger RNA
        self.trigRNA = [[[~dC[j], ~dB[j], ~dA[j]] + term2[j] for j in range(self.N)]]

    def reverse_toehold_switch(self, dimensions=[10,10,10,10,10,6]):
        '''
        Sets design specifications for ts45 (reverse toehold switch)
        OFF -> ON trigger logic
        g1 + t1 -> t1g1
        '''
        # unpack dimensions
        [A,B,C,D,E,F,G,L] = dimensions

        # fixed domains
        target = self.target
        DU_cas9h = self.DU_cas9h
        cas9h = self.cas9h 
        term2 = self.term2
        g1Tm = self.g1Tm
        g2Tm = self.g2Tm
        t1Tm = self.t1Tm
        DU_g1Tm = self.DU_g1Tm
        DU_g2Tm = self.DU_g2Tm
        DU_t1Tm = self.DU_t1Tm
        
        # variable domains
        dA = [nupack.Domain('N'*A, name=['A',j]) for j in range(self.N)]
        dB = [nupack.Domain('N'*B, name=['B',j]) for j in range(self.N)]
        dC = [nupack.Domain('N'*C, name=['C',j]) for j in range(self.N)]
        dD = [nupack.Domain('N'*D, name=['D',j]) for j in range(self.N)]
        dE = [nupack.Domain('N'*E, name=['E',j]) for j in range(self.N)]
        dF = [nupack.Domain('N'*F, name=['F',j]) for j in range(self.N)]
        dG = [nupack.Domain('N'*G, name=['G',j]) for j in range(self.N)]
        dL = [nupack.Domain('N'*L, name=['L',j]) for j in range(self.N)]
        # specify strands
        term2 = [[dE[j], dF[j], ~dE[j], dG[j]] for j in range(self.N)]
        DU_term2 = ['D'+str(dE[j].nt())+'(U'+str(dF[j].nt())+')U'+str(dG[j].nt()) for j in range(self.N)]

        Sg1 = [nupack.TargetStrand([dB[j], target[j], cas9h, dL[j], dC[j], dD[j], ~dC[j]] + term2[j] + [~dB[j], ~dA[j], g1Tm], name=['spec_g1',j]) for j in range(self.N)]
        St1 = [nupack.TargetStrand([dA[j], dB[j], t1Tm], name=['t1',j]) for j in range(self.N)]
 
        Cg1 = [nupack.TargetComplex([Sg1[j]], structure='D'+str(dB[j].nt())+'(U'+str(target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+'D'+str(dC[j].nt())+'(U'+str(dD[j].nt())+')'+DU_term2[j]+')U'+str(dA[j].nt())+DU_g1Tm,
                name=['Cg1', j]) for j in range(self.N)]
        Ct1 = [nupack.TargetComplex([St1[j]],
                structure='U'+str(dA[j].nt()+dB[j].nt())+DU_t1Tm,
                name=['Ct1',j]) for j in range(self.N)]
        Ct1g1 = [nupack.TargetComplex([Sg1[j], St1[j]],
                structure='U'+str(dB[j].nt()+target[j].nt())+DU_cas9h+'U'+str(dL[j].nt())+'D'+str(dC[j].nt())+'(U'+str(dD[j].nt())+')'+DU_term2[j]+'D'+str(dB[j].nt()+dA[j].nt())+'('+DU_g1Tm+'+)'+DU_t1Tm,
                name=['Ct1g1', j]) for j in range(self.N)]
        
        # test tubes
        steps = []
        for j in range(self.N):
            # g1 -> g1
            steps.append(nupack.TargetTube({Cg1[j]:1e-8}, nupack.SetSpec(max_size=2), name=['step0',j]))
            # g1 + t1 -> t1g1
            steps.append(nupack.TargetTube({Ct1[j]:1e-8, Cg1[j]:1e-8}, nupack.SetSpec(max_size=2, exclude=[Ct1g1[j]]), name=['step1',j]))
            # t1g1 -> t1g1
            steps.append(nupack.TargetTube({Ct1g1[j]:1e-8}, nupack.SetSpec(max_size=2), name=['step2',j]))
        self.tubes = steps
        
        # set the weights
        weights = nupack.Weights(self.tubes)
        for j in range(self.N):
            weights[target[j], :, :, :] *= 0
        self.weights = weights

        # define domain scopes for soft and hard constraints
        self.scope = [[dL[j], dC[j], dD[j], ~dC[j]] + term2[j] + [~dB[j], ~dA[j]] for j in range(self.N)]
        
        # define trigRNA domains
        self.trigRNA = [[[dA[j], dB[j]] for j in range(self.N)]]

