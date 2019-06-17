# Fraction of Reads in Peaks (FRiP)
Calculating FRiP for ChIP-seq or ChEC-seq data at various downsampled read numbers.

# About
This workflow uses software from a containerized conda environment to downsample reads, call peaks, and then get number of reads in peaks. This information is useful in assessing background signal, as well as sequencing saturation for ChIP-seq and ChEC-seq samples.

# Getting Started

## Cloning Repository

To get started, you must first clone the FRiP repository. Navigate to a directory you would like to clone the repository to and enter `git clone https://github.com/rpolicastro/chip_downsampling.git`.

## Installing Singularity

Singularity containers are self contained 'boxes' that house the software and other files necessary for the workflow. The container itself will automatically be downloaded, but you must have the Singularity software installed to both download and use the container. Please refer to the [documentation](https://www.sylabs.io/docs/) on their website.

## Specifying Run Settings

The last step is to set a few settings in the 'settings.conf' file in the main repository directory. An example settings file is provided in the 'DOCS' directory of the repository.

 Setting | Description |
| ------- | ----------- |
|REPDIR|Main repository directory.|
|WORDKIR|Directory where outputs will be saved.|
|BAM|Path and file name of input BAM.|
|PAIRED|Whether the BAM is paired end or not (TRUE/FALSE).|
|FROM|Starting number of reads to downsample to.|
|TO|Ending number of reads to downsample to.|
|BY|Intervals to downsamples to between the FROM and TO values.|
|CONTROL|Input/control BAM file path and file name.|
|GENOME_SIZE|Effective genome size.|

## Running the Workflow

After getting singularity installed, and the settings specified, you are now ready to run the workflow. Navigate to the main directory and enter 'bash main.sh'.

##### Notes for IU Folks
If you wish to submit the workflow to a compute node, you can do so by submitting it through the TORQUE resource manager. Navigate to the directory that contains both your 'settings.conf' and 'main.sh' files. Create a file called **submit_workflow.sh** with the following contents:

```
#!/bin/bash

## Navigate back to directory containing the 'main.sh' and 'settings.conf' file.
cd $PBS_O_WORKDIR

## Load the singularity module on Carbonate/Karst
module load singularity/3.2.0

## Start the workflow
bash main.sh
```

You can now submit the workflow.`qsub -l nodes=1:ppn=8,vmem=64gb,walltime=12:00:00 submit_workflow.sh`. 'ppn' specifies the threads/cores, and 'vmem' is the virtual memory.

# Built With

This workflow would not be possible without the great software listed below.

- [Anaconda](https://www.anaconda.com/) - Software package manager and virtual environment.
- [Singularity](https://www.sylabs.io/docs/) - Containerize sofware and files.
- [Samtools](http://www.htslib.org/) - SAM/BAM manipulation.
- [MACS](https://github.com/taoliu/MACS) - Peak caller.
- [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) - Annotating fragments to overlapping peaks.
