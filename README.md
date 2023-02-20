# CGDesigner

`CGDesigner` is a python based tool for designing and analyzing conditional guide RNAs with NUPACK. This library and python executable is organized into the following submodules:

`design` contains the commands to generate nupack design specifications which are used as input in the checkpoint module

`checkpoint` contains the commands to run the nupack designs locally and generate cgRNA sequences

`ftp` contains the commands to transfer files to and from the s3 file server

`unpack` unpacks nupack design output and data fields into csv format

`kube` contains the commands to send jobs to the kubernetes computing cluster

`analysis` contains the commands to filter and perform test tube analysis on cgRNA sequences

`oligos` contains commands to generate oligos for cloning

# Installation
This package requires `nupack` to be installed on the server. This package can be found at [nupack.org](https://nupack.org/download/overview).

Run the following commands to clone and install this tool.
```
git clone https://github.com/zchen15/CGDesigner`
cd CGDesigner
pip3 install .
```

# Usage 



# References
If you use `CGDesigner` in a publication, please cite:

Hanewich-Hollatz MH, Chen Z, Hochrein LM, Huang J, Pierce NA. Conditional Guide RNAs: Programmable Conditional Regulation of CRISPR/Cas Function in Bacterial and Mammalian Cells via Dynamic RNA Nanotechnology. ACS Cent Sci. 2019 Jul 24;5(7):1241-1249. doi: [10.1021/acscentsci.9b00340](https://doi.org/10.1021/acscentsci.9b00340).

