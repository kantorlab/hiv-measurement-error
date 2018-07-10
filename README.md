Scripts to reproduce results in the manuscript:

M Howison, M Coetzer, R Kantor. 2018. Measurement error and variant-calling in
deep Illumina sequencing of HIV.  *bioRxiv* 276576.
doi:[10.1101/276576](https://doi.org/10.1101/276576)

## Install

An Ubuntu 16.04-based Docker image with all prequisites installed is available
from [DockerHub](https://hub.docker.com/r/kantorlab/hiv-measurement-error). *If
you choose to use the Docker image, please read the notes about resource
requirements in the "Running" section below. In particular, you may need to
adjust the memory settings in Docker to allow for at least
10GB of memory.*

The scripts have been tested on CentOS 6.8 and Ubuntu 16.04, but may work on
other distributions of 64-bit Linux.

### Conda packages

Compiled Linux 64-bit binaries are available for all dependencies in the
[kantorlab Anaconda channel](https://anaconda.org/kantorlab). The corresponding
conda recipes are available from
[https://github.com/kantorlab/conda-recipes](https://github.com/kantorlab/conda-recipes).

To use the packages, first install [Anaconda 3](https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh).
Then recreate the analysis environment with the command:

    conda create -n hiv-measurement-error -c kantorlab blastn=2.7.1 hivmmer=0.1.2 iva=1.0.9 mafft=7.313 matplotlib quasitools=0.3.1 sra-tools=2.9.1.1 scons=3.0.1.1 trimmomatic=0.36
    conda install -n hiv-measurement-error seaborn

To activate the environment, use:

    source activate hiv-measurement-error

### Manual installation

Alternatively, the following list of prerequisites can be installed manually.
**Note**: scons was patched to support python 3.6.5 and is available in a fork
[https://github.com/mhowison/scons/releases](https://github.com/mhowison/scons/releases)
as release "3.0.1-hotfix1".

* [bambamc](https://github.com/gt1/bambamc) 0.0.50
* [biopython](https://biopython.org/) 1.71
* [blast](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs) 2.7.1
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.3.4.1
* [fastaq](https://github.com/sanger-pathogens/Fastaq) 3.15.0
* [fastx_toolkit](https://github.com/agordon/fastx_toolkit) 0.0.14
* [hivmmer](https://github.com/kantorlab/hivmmer) 0.1.2
* [hmmer](http://hmmer.org/) 3.2.1
* [iva](https://github.com/sanger-pathogens/iva) 1.0.9
* [kmc](https://github.com/refresh-bio/KMC) 3.0.0
* [mafft](http://mafft.cbrc.jp/alignment/software/) 7.313
* [matplotlib](https://matplotlib.org/) 2.2.2
* [mummer](http://mummer.sourceforge.net) 3.23
* [numpy](http://www.numpy.org/) 1.14.5
* [pandas](https://pandas.pydata.org/) 0.23.1
* [pear](http://www.exelixis-lab.org/web/software/pear) 0.9.11
* [python](https://www.python.org) 3.6.5
* [pysam](https://github.com/pysam-developers/pysam) 0.14.1
* [quasitools](https://github.com/phac-nml/quasitools) 0.3.1 (contains HyDRA)
* [samtools](https://github.com/samtools/samtools) 1.3.1
* [scipy](https://www.scipy.org/) 1.1.0
* [scons](https://github.com/mhowison/scons/releases) 3.0.1-hotfix1
* [seaborn](https://seaborn.pydata.org) 0.8.1
* [smalt](http://www.sanger.ac.uk/science/tools/smalt-0) 0.7.6
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 0.36

Additionally, the following standard Unix utilities and OS packages are needed
(Ubuntu package names are specified):
* bc
* bzip2
* git
* openjdk-9-jre-headless (or an equivalent Java JRE)
* unzip
* util-linux (for the `getopt` command)
* wget

## Running

The run order and dependencies of the scripts are specified in the SConstruct
file.  The entire analysis can be run by executing the `scons` command from the
root directory of the repo.

Several steps require extensive compute and memory resources, and the analysis is best
run on a compute cluster, where it can be parallelized by SCons. To run the analysis
on a compute cluster that uses the [SLURM](https://slurm.schedmd.com/) batch
system, you can simply uncomment the following line in SConstruct:

    #srun = "srun"

and then execute the analysis with multiple job slots with `scons -j32` (e.g. for
scheduling up to 32 jobs at a time in SLURM).  This will prefix many of the
commands with the `srun` command so that they are scheduled and launched
through the SLURM batch system.  If your cluster requires more specific
arguments to `srun` (such as a partition or account name), you can append these
to the `srun` variable in SConstruct.

Overall, the entire analysis requires approximately 100 CPU-hours and the
maximum memory required for any step is approximately 10GB. **On a single
workstation, this means the entire analysis could take around a week to run**.
Therefore, we recommend you run the analysis on a compute cluster if possible.
The majority of compute time is required for pidalign/pidalyse and IVA/shiver,
and you may see long pauses during those steps.

### Scratch directory

The analysis will generate many large intermediate files that should be written
to a file system with good performance and capacity. All of these files are
placed in the "scratch" subdirectory, which can be a symlink to another file system.
On a compute cluster, you may want to clone the repo to your home directory,
but locate the scratch directory on another high-performance file system. For
example, on our compute cluster, we achieve this with:

    cd $HOME
    git clone https://github.com/kantorlab/hiv-measurement-error
    cd hiv-measurement-error
    mkdir -p /gpfs/scratch/mhowison/hiv-measurement-error
    ln -s /gpfs/scratch/mhowison/hiv-measurement-error scratch

where `/gpfs/scratch` is a high-performance file system.

The total scratch storage used by the entire analysis is around 40GB.

### Results

A copy of the results are included in the git repository in the `results`
subdirectory.

Instead of running the entire analysis, an individual result can be recomputed
by calling scons on a specific target.  For example, to recreate the coverage
plot for 5VM, use:

    scons results/coverage.eps

This will automatically build any missing intermediate files in the scratch
directory.
