# Fraction of Reads in Peaks (FRiP)
Calculating FRiP for ChIP-seq or ChEC-seq data at various downsampled read numbers.

# About
This workflow uses software from a conda environment to downsample reads, call peaks, and then get number of reads in peaks. This information is useful in assessing background signal, as well as sequencing saturation for ChIP-seq and ChEC-seq samples.

# Getting Started

## Cloning Repository

To get started, you must first clone the FRiP repository. Navigate to a directory you would like to clone the repository to and enter `git clone https://github.com/rpolicastro/chip_downsampling.git`.

## Preparing Conda Environment

This workflow takes advantage of the [conda](https://conda.io/en/latest/) package manager and virtual environment. The conda package manager installs both the main software and all dependencies into a 'virtual environment' to ensure compatabilty. Additionally, the provided 'environments.yml' file can be used to install the same major software versions as used to develop the script, ensuring prolonged compatability and reproducibility.

Before creating the environment, you must first install miniconda.
1. [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda), and make sure that conda is in your PATH.
2. Update conda to the latest version `conda update -n base -c defaults conda`.

You are now ready to create the virtual sofware environment, and download all software and dependencies. It is recommended to reproduce the environment used when creating the script, but instructions on installing the latest software are provided as an alternative.

#### Reproducing the Development Environment (Recommended)

To install the major software versions used when developing this script, navigate to the 'DOCS' directory, and use the provided 'environments.yml' file to create your conda environment.
```
conda env create -f environment.yml
```
For posterity, *all* software and versions used when developing the script are provided in the 'development_environment.yml' file located in the 'DOCS' directory for the repository. This file can *not* be used to install the environment on your computer, because many of the dependencies and software builds are system specific. However, this file may help you troubleshoot any dependency errors that may occur in your environment.

#### Installing The Latest Software Versions

1. Create the new environment and specify the software to include in it.
```
conda create -n chip-downsampling -c conda-forge -c bioconda \
r-getopt samtools bioconductor-rsubread macs2
```
2. Update the software to the latest compatible versions.
```
conda update -n chipseq-automation -y -c conda-forge -c bioconda --all
```

If you wish to use any of the software in the environment outside of the workflow you can type `conda activate chip-downsampling`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.

## Specifying Run Settings

The last step is to set a few settings in the 'settings.conf' file in the main repository directory. An example settings file is provided in the 'DOCS' directory of the repository.

 Setting | Description |
| ------- | ----------- |
|WORDKIR|Directory where outputs will be saved.|
|BAM|Path and file name of input BAM.|
|PAIRED|Whether the BAM is paired end or not (TRUE/FALSE).|
|FROM|Starting number of reads to downsample to.|
|TO|Ending number of reads to downsample to.|
|BY|Intervals to downsamples to between the FROM and TO values.|
|CONTROL|Input/control BAM file path and file name.|
|GENOME_SIZE|Effective genome size.|

## Running the Workflow

After getting the conda environment ready, the sample sheet prepared, and the settings specified, you are now ready to run the workflow. Navigate to the main directory and enter 'bash main.sh'.

# Built With

This workflow would not be possible without the great software listed below.

- [Anaconda](https://www.anaconda.com/) - Software package manager and virtual environment.
- [Samtools](http://www.htslib.org/) - SAM/BAM manipulation.
- [MACS](https://github.com/taoliu/MACS) - Peak caller.
- [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) - Annotating fragments to overlapping peaks.
