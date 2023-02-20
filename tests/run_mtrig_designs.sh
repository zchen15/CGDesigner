#!/usr/bin/env bash
CGDesigner -c npbop.json design -material rna -s rna_trig -N 4 -o rnatrig -d 20 -fstop 0.01

CGDesigner -c npbop.json design -material rna -s reverse_toehold_switch -N 1 -gin ../data/mtrig.csv -o ts45_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 50 10
CGDesigner -c npbop.json design -material rna -s terminator_switch -N 1 -gin ../data/mtrig.csv -o ts32_mtrig -d 10 20 6 0 7 3 1 -fstop 0.01 -scan 100 100
CGDesigner -c npbop.json design -material rna -s split_guideRNA -N 1 -gin ../data/mtrig.csv -o ts26_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 100 100
CGDesigner -c npbop.json design -material rna -s inhibited_split_terminator_switch -N 1 -gin ../data/mtrig.csv -o ts38_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 100 100
CGDesigner -c npbop.json design -material rna -s reverse_split_terminator_switch -N 1 -gin ../data/mtrig.csv -o ts25_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 100 100
CGDesigner -c npbop.json design -material rna -s single_hairpin_switch -N 1 -gin ../data/mtrig.csv -o ts18b_mtrig -d 10 20 10 7 3 1 6 -fstop 0.01 -scan 100 100
CGDesigner -c npbop.json design -material rna -s two_hairpin_switch -N 1 -gin ../data/mtrig.csv -o ts23b_mtrig -d 10 20 10 0 0 0 6 0 7 3 1 -fstop 0.01 -scan 100 100
CGDesigner -c npbop.json design -material rna -s catalytic_two_hairpin_switch -N 1 -gin ../data/mtrig.csv -o ts12_mtrig -d 10 20 10 6 -fstop 0.01 -scan 100 100

