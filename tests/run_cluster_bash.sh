CGDesigner -c data/npbop.json design -material rna -s reverse_toehold_switch -N 1 -gin data/mtrig.csv -o ts45_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 100 100
CGDesigner -c data/npbop.json checkpoint -i ts45*.spec

