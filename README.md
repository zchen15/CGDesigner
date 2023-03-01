# CGDesigner
`CGDesigner` is a python based tool for designing and analyzing conditional guide RNAs with NUPACK. This library and python executable is organized into the following submodules:

`design` contains the commands to generate nupack design specifications which are used as input in the checkpoint module

`checkpoint` contains the commands to run the nupack designs locally and generate cgRNA sequences

`ftp` contains the commands to transfer files to and from the s3 file server

`unpack` unpacks nupack design output and data fields into csv format

`kube` contains the commands to send jobs to the kubernetes computing cluster

`analysis` contains the commands to filter and perform test tube analysis on cgRNA sequences

`oligos` contains commands to generate oligos for cloning

## Installation
This package requires `nupack` to be installed on the server. This package can be found at [nupack.org](https://nupack.org/download/overview).

Run the following commands to clone and install this tool.
```
git clone https://github.com/zchen15/CGDesigner
cd CGDesigner
pip3 install .
```

## Usage 
The `tests/` directory contains many example scripts for generating conditional guide RNAs. The following are a few useful excerpts to help you start off.

### Generating NUPACK designs
`design` is used to generate nupack design specifications which can be run locally or on the kubernetes cluster. The following are a few examples designs specifications that can be generated.

```
# Get help on design options
CGDesigner design -h

# Generate 4 orthogonal universal trigger RNAs of length 20nt
CGDesigner design -material rna -s rna_trig -N 4 -o rnatrig -d 20 -fstop 0.01

# Generate reverse toehold switch cgRNAs targeting trigger sequences in mtrig.csv
CGDesigner design -material rna -s reverse_toehold_switch -gin data/mtrig.csv -o ts45_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 50 10

# Generate terminator switch cgRNAs targeting trigger sequences in mtrig.csv
CGDesigner design -material rna -s terminator_switch -gin data/mtrig.csv -o ts32_mtrig -d 10 20 6 0 7 3 1 -fstop 0.01 -scan 100 100
```
`-o` defines the output file prefix

`-scan` generates designs targeting 50nt sequences at stride 10nt for the sequences in mtrig.csv

`-d` defines the domain dimensions to use which vary from design to design

`-gin` defines the trigger RNA input to design against, if none are provided then trigger sequenes are considered unconstrained

`-fstop` defines the threshold ensemble defect to reach for design to stop

`-material` defines the thermodynamic model to use

`-s` defines the test_tube formula to use. These are functions listed in `formula.py`

Details about these designs can be found in [my thesis](https://github.com/zchen15/CaltechThesis/raw/main/revisions/Thesis_20230225ZC.pdf). New design formulas can be added to `CGDesigner/formula.py` and called natively from CGDesigner via `-s`.

### Running jobs locally
Design jobs can be run locally by providing design `.spec` files as input. These are equivalent to design checkpoints that can be restarted anytime a job fails.

```
# Get help on this submodule
CGDesigner checkpoint -h

# Runs designs with trials = 4 and checkpoint interval = 30
CGDesigner -c data/npbop.json checkpoint -i *.spec -trials 4 -cint 30 
```

`-trials` defines number of jobs to running in parallel

`-cint` defines the checkpoint interval in seconds

`-c` defines the path to the config file which has other default run parameters

### Submitting jobs to the cluster
The following example shows how design jobs can be submitted to the kubernetes cluster using the `kube` submodule. Make sure you are using the same version of NUPACK as that installed on the cluster otherwise the design `.spec` files will fail to load. This submodule can also be used to submit and run bash scripts to the docker container.

```
# get help on this submodule
CGDesigner kube -h

# submits design jobs to the cluster
CGDesigner -c data/npbop.json kube -i *.spec

# submits bash script to run on the cluster
CGDesigner -c data/npbop.json kube -i test.sh

# lists current running jobs
CGDesigner -c data/npbop.json kube -ls '*'

# lists current running pods
CGDesigner -c data/npbop.json kube -ls '*' -pods

# removes jobs matching keyword ts45*
CGDesigner -c data/npbop.json kube -rm 'ts45*'

# removes pods matching keyword ts45*
CGDesigner -c data/npbop.json kube -i *.spec -pods

# clears completed or failed pods and jobs
CGDesigner -c data/npbop.json kube -clear
CGDesigner -c data/npbop.json kube -clear -pods
```

`-c` defines the path to the config file. This files contains information about s3 and kubernetes cluster credentials.

### Download finished results from the s3 file server
The following shows how to download results from the s3 file server using the `ftp` submodule.
```
# get help on this submodule
CGDesigner ftp -h

# get list of ts45 results
CGDesigner -c data/npbop.json ftp -ls 'ts45*.csv'

# Download ts45 results
CGDesigner -c data/npbop.json ftp -get 'ts45*.csv'

# Remove ts45 results
CGDesigner -c data/npbop.json ftp -rm 'ts45*.csv'
```

### Filtering designs
The following example shows how to filter and analyze cgRNA sequences with the `analysis` submodule.
```
# download results and filter by prediction
CGDesigner -c ../data/npbop.json ftp -get 'ts45_mRNA_1*.csv'
# note -noterm strips the terminator sequence
CGDesigner analysis -i ts45_mRNA_0*.csv -o strands.csv -m remove_homopolymers  -noterm

# select designs targeting the last 400 bp of the PAX7 gene
# -s selects for Strands with keyword '*t1*'
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
```

### Oligo generation
The following shows how to generate oligos with BsaI golden gate sites padded to 300nt. The output can be upload to Twist or IDT to obtain gene block or oligo pool order POs.
```
# output is written to oligos.csv
# -g is used to add golden gate sites or flanking sequences
# -pad defining the padding length. Here all strands are made to be at least 300nt
echo 'generating oligos for the designs'
CGDesigner oligos -noterm -m twist_outer -g "tatatagGGTCTCcCACA " " CTTTgGAGACCctatata" -s "*g1*0*" -pad 300 -i strands.csv -o oligos.csv
```

### Rebuilding the docker container
The docker contain can be rebuild and added to the kubernetes cluster by running `bash docker_build.sh`. This will update the container with new code you added to this directory.


## Issues
If you experience any issues with the code, please post them on the issues section along with the log file. I will monitor this periodically and try to fix issues as they arise.

## References
If you use `CGDesigner` in a publication, please cite:

Hanewich-Hollatz MH, Chen Z, Hochrein LM, Huang J, Pierce NA. Conditional Guide RNAs: Programmable Conditional Regulation of CRISPR/Cas Function in Bacterial and Mammalian Cells via Dynamic RNA Nanotechnology. ACS Cent Sci. 2019 Jul 24;5(7):1241-1249. doi: [10.1021/acscentsci.9b00340](https://doi.org/10.1021/acscentsci.9b00340).

