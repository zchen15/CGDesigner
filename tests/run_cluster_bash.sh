CGDesigner -c data/npbop.json design -material rna -s terminator_switch -N 1 -gin data/mtrig.csv -o ts32_mtrig -d 10 20 6 0 7 3 1 -fstop 0.01 -scan 100 100
CGDesigner -c data/npbop.json checkpoint -i ts32*.spec

