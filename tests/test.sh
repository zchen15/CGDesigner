#!/bin/bash

# A collection of bash scripts to test the functionality of the CGDesigner tool

# test generation of design spec
#bash run_mtrig_designs.sh

# test checkpoint and design start locally
if [ $1 == checkpoint ];then
    CGDesigner -c ../data/npbop.json checkpoint -i *.spec -trials 4 -cint 30 
fi

if [ $1 == server ];then
    # test file transfer
    CGDesigner -c ../data/npbop.json ftp -ls '*'

    # test sending designs to kubernetes 
    CGDesigner -c ../data/npbop.json kube -i *.spec
CGDesigner -c ../data/npbop.json kube -ls '*'
fi

# download results and filter by prediction
CGDesigner -c ../data/npbop.json ftp -get 'ts45_mRNA_1*.csv'
CGDesigner analysis -i ts45_mRNA_0*.csv -o strands.csv -m remove_homopolymers  -noterm

# select designs targeting the last 400 bp of the PAX7 gene
CGDesigner analysis -i strands.csv -o strands.csv -m add_position -r ../data/PAX7_400.csv -s '*t1*'
CGDesigner analysis -i strands.csv -o strands.csv -m add_position -r ../data/PAX7.csv -s '*t1*'

# filter by test tube specifications to make sure they still work
echo 'running test tube analysis on designs'
CGDesigner -v analysis -material rna -i strands.csv -o strands.csv -m filter_ts45 -r ../data/PAX7.csv

# filter for best 100 designs based on external prediction score
echo 'filtering by external predictor'
CGDesigner analysis -i strands.csv -o strands.csv -m add_prediction -r ../data/prediction.csv
CGDesigner analysis -i strands.csv -o strands.csv -m get_n_best -d 100 -30 -r 'prediction'

# filter for best 100 designs based on design defect
echo 'filtering by design defect'
CGDesigner analysis -i strands.csv -o strands.csv -m parse_defect
CGDesigner analysis -i strands.csv -o strands.csv -m get_n_best -d 100 -20 -r 'defect'

# get the 20 best designs based on external prediction score
echo 'filter for best 20 designs based on external predictor'
CGDesigner analysis -material rna -i strands.csv -o strands.csv -m get_n_best -d 16 -20 -r 'prediction'


# generate bar plot of off and on target mRNAs



# generate oligos from filtered designs
echo 'generating oligos for the designs'
CGDesigner oligos -noterm -m twist_outer -g "tatatagGGTCTCcCACA " " CTTTgGAGACCctatata" -s "*g1*0*" -pad 300 -i strands.csv -o oligos.csv

# generate plasmid maps


