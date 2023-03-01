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
The `tests/` directorry contains many example scripts for generating conditional guide RNA. The following are a few useful excerpts.

### Generating NUPACK designs
`design` is used to generate nupack design specifications which can be run locally or on the kubernetes cluster. The following are a few examples designs specifications that can be generated.

```
# Get help on design options
CGDesigner design -h

# Generate N orthogonal universal trigger RNAs
CGDesigner design -material rna -s rna_trig -N 4 -o rnatrig -d 20 -fstop 0.01

# Generate reverse toehold switch cgRNAs targeting trigger sequences in mtrig.csv
CGDesigner design -material rna -s reverse_toehold_switch -gin data/mtrig.csv -o ts45_mtrig -d 10 20 10 3 10 3 3 6 -fstop 0.01 -scan 50 10

# Generate terminator switch cgRNAs targeting trigger sequences in mtrig.csv
CGDesigner design -material rna -s terminator_switch -gin ../data/mtrig.csv -o ts32_mtrig -d 10 20 6 0 7 3 1 -fstop 0.01 -scan 100 100
```

`-scan` generates designs targeting 50nt sequences at stride 10nt for the sequences in mtrig.csv
`-d` defines the domain dimensions to use which vary from design to design
`-gin` defines the trigger RNA input to design against, if none are provided then trigger sequenes are considered unconstrained
`-fstop` defines the threshold ensemble defect to reach for design to stop
`-material` defines the thermodynamic model to use

Details about these designs can be found in [my thesis](https://github.com/zchen15/CaltechThesis/raw/main/revisions/Thesis_20230225ZC.pdf). New design formulas can be added to `CGDesigner/formula.py`

### Running jobs locally
Design jobs can be run locally by providing design `.spec` files as input. These are equivalent to design checkpoints that can be restarted anytime a job fails.

```
CGDesigner checkpoint -i *.spec -trials 4 -cint 30 
```

`-trials` defines number of jobs to running in parallel
`-cint` defines the checkpoint interval in seconds

### Submitting jobs to the cluster
Design jobs can be submitted to the kubernetes cluster using the following. Make sure you are using the same version of NUPACK as that installed on the cluster otherwise the design `.spec` files will fail to load.

```
# submits design jobs to the cluster
CGDesigner -c ../data/npbop.json kube -i *.spec

# lists current running jobs
CGDesigner -c ../data/npbop.json kube -ls '*'

# lists current running pods
CGDesigner -c ../data/npbop.json kube -ls '*' -pods

# removes jobs matching keyword ts45*
CGDesigner -c ../data/npbop.json kube -rm 'ts45*'

# removes pods matching keyword ts45*
CGDesigner -c ../data/npbop.json kube -i *.spec -pods

# clears completed or failed pods and jobs
CGDesigner -c ../data/npbop.json kube -clear
CGDesigner -c ../data/npbop.json kube -clear -pods
```

`-c` defines the path to the config file. This files contains information about s3 and kubernetes cluster credentials.

### Filtering designs


### Rebuilding the docker container
The docker contain can be rebuild and added to the kubernetes cluster by running `bash docker_build.sh`. This will update the container with new code you added to this directory.

## References
If you use `CGDesigner` in a publication, please cite:

Hanewich-Hollatz MH, Chen Z, Hochrein LM, Huang J, Pierce NA. Conditional Guide RNAs: Programmable Conditional Regulation of CRISPR/Cas Function in Bacterial and Mammalian Cells via Dynamic RNA Nanotechnology. ACS Cent Sci. 2019 Jul 24;5(7):1241-1249. doi: [10.1021/acscentsci.9b00340](https://doi.org/10.1021/acscentsci.9b00340).

