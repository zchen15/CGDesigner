import subprocess

# sort through the five ts45 design dimensions
for i in range(0,6):
    cmd = []
    cmd+= ['CGDesigner -c ../data/npbop.json ftp -get "ts45_mRNA_'+str(i)+'*.csv" ']
    cmd+= ['CGDesigner analysis -i ts45_mRNA_'+str(i)+'*.csv -o strands.csv -m remove_homopolymers -noterm -s "*g1*" -ex ggtctc gagacc']
    #cmd+= ['CGDesigner analysis -i strands.csv -o strands.csv -m add_position -r ../data/PAX7_400.csv -s "*t1*"']
    cmd+= ['CGDesigner analysis -i strands.csv -o strands.csv -m add_position -r ../data/PAX7.csv -s "*t1*"']
    cmd+= ['CGDesigner -v analysis -material rna -i strands.csv -o strands.csv -m filter_ts45 -r ../data/PAX7.csv']
    cmd+= ['CGDesigner analysis -i strands.csv -o strands.csv -m add_prediction -r ../data/prediction.csv']
    cmd+= ['CGDesigner analysis -i strands.csv -o strands.csv -m get_n_best -d 100 -30 -r "prediction"']
    cmd+= ['CGDesigner analysis -i strands.csv -o strands.csv -m parse_defect']
    cmd+= ['CGDesigner analysis -i strands.csv -o strands.csv -m get_n_best -d 100 -20 -r "defect"']
    cmd+= ['CGDesigner analysis -material rna -i strands.csv -o strands_'+str(i)+'.csv -m get_n_best -d 16 -20 -r "pos1"']
    for c in cmd:
        subprocess.call(c, shell=True)

        