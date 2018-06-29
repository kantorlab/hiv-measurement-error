Scripts to reproduce results from the manuscript:

*Measurement error and variant-calling in deep Illumina sequencing of HIV*

## Conda environment

Compiled Linux 64-bit binaries are available for all dependencies in the
[kantorlab Anaconda channel](https://anaconda.org/kantorlab).

To use the packages, first install [Anaconda 3](https://conda.io/miniconda.html).
Then recreate the analysis environment with the command:

    conda create -n hiv-measurement-error -c kantorlab blastn=2.7.1 hivmmer=0.1.1 iva=1.0.9 mafft=7.313 matplotlib quasitools=0.3.1 scons=3.0.1.1 trimmomatic=0.36

To activate the environment, use:

    source activate hiv-measurement-error

## Scratch directory

The analysis will generate many large intermediate files that should be written
to a file system with good performance and capacity. All of these files are
placed in the "scratch" subdirectory, which you need to create. On our HPC
system we create it as a link to a directory in a scratch file system like
this:

    mkdir -p /gpfs/scratch/mhowison/hiv-measurement-error
    ln -s /gpfs/scratch/mhowison/hiv-measurement-error scratch

## Results

A copy of the results are included in the git repository in the `results`
subdirectory.

Each result can be recomputed by calling biomake (included in the Anaconda
channel) on the appropriate target in the included Makefile. For example, to
recreate
the coverage plot for 5VM, use:

    biomake results/coverage.eps

This will automatically build any missing intermediate files in the scratch
directory.
