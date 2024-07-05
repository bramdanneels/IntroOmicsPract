# Practical 0: Introduction and Setup

This practical contains general information necessary to complete the practicals, and information on how to set up certain software necessary for the practicals.

## IGV
In practical 1 and onwards, you will be using a genome browser to visualise genome data. We will be using [IGV] for this. IGV is an easy-to-use, java-based genome browser. It can automatically download a range of reference genomes and their annotations, but can also be used for visualisation of your own data.
IGV can be downloaded from [this link](https://igv.org/doc/desktop/#DownloadPage/)

## NREC & Linux
Praticals 2-6 will be using a range of -omics tools. Most -omics/bioinformatic tools are developped for linux, and require some computational resources. 
Students of the BINF201 course will have access to an [NREC](https://www.nrec.no/) instance with the necessary computational power and pre-installed software. Please check the BINF201 MittUiB course page for more information on how to get access to the BIN201 server.
It is also possible to do the practicals on your own computer. If you have a Linux system, you can just install the necessary software and download the necessary data yourself (instructions will be provided in the tutorials themselves).
If you have a Windows system, it is advised to use [WSL (Windows subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install). For Mac users, it can be possible to run most of the software on the command line as well, but this has not been tested (and has caused problems in previous iterations of the course).

## Shell
Most software will be controlled from the command line (the shell). If you are not familiar with working on the command line, please have a look at [this excellent shell tutorial](https://swcarpentry.github.io/shell-novice/).
A basic mastery of the command line will make going through the tutorials and identifying possible problems a lot easier.

## Conda
The necessary software is already installed on the NREC server. Most software have been installed using the [Conda](https://docs.conda.io/en/latest/) package manager. Conda allows users to easily install a range of (bio)informatics software. The main benefit of using Conda is that it automatically tries to solve dependencies.
This means that if you ask Conda to install a certain software package (or python package), it will recursively check which other packages are required, and install them. 
Software in this course has been installed using [Miniforge](https://github.com/conda-forge/miniforge). Miniforge has two main benefits over using the standard conda:
- It uses the faster and more powerfull [Mamba](https://github.com/mamba-org/mamba) package manager, which is basically an upgrade of Conda. 
- It preconfigures conda-forge as a package repository, allowing default access to a whole range of extra packages.
Many of the bioinformatic packages that we will use, are only available through [Bioconda](https://bioconda.github.io/index.html). If you are installing the packages yourself, it is advised to set up the bioconda channel as part of the default by running the following code:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## Data download

All data used in these practicals is publically available. The data is provided on the NREC server for the students of the course. 

For people who do not have access to NREC and/or want to do the practicals on their own computer, all data used in this practicals is available through [Zenodo]().
Download links to the specific files needed for a tutorial are available in each tutorial.

Alternatively, you can download and pre-process the data yourself. To do this, it is advised to create an environment dedicated to downloading and processing the necessary data.
You can create this enviroment by running the following command (after having installed `mamba`, see above):

```
mamba create -n downloader sra-tools ncbi-datasets-cli seqtk pigz
```

- [SRA tools](https://github.com/ncbi/sra-tools) is used for downloading data from [SRA](https://www.ncbi.nlm.nih.gov/sra) (mainly raw sequencing data in `.fastq` format.
- [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) is used for downloading assembled sequences (genomes, genes, proteins) from NCBI.
- [SeqTK](https://github.com/lh3/seqtk) is a sequence toolkit, which we will mainly use to subset larger datasets into smaller ones
- [Pigz](https://github.com/madler/pigz) is a parallel implementation of GZip, allowing fast (de)compression of files

The commands needed to download the necessary data are included in each practical as well.

## R and Rstudio
For practicals 7-9 we will be using R and R studio. These have to be downloaded and installed yourself. Please follow the instructions [here](https://posit.co/download/rstudio-desktop/) on how to install both R and Rstudio.
For these practicals, it is not necessary to be on Linux or the NREC server.
