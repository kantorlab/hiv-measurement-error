#!/bin/bash
set -e
cd scratch
prefetch $1
fastq-dump --defline-qual "+" --split-files --defline-seq '@$sn[_$rn]/$ri' $1
mv $1_1.fastq $2_1.fastq
mv $1_2.fastq $2_2.fastq
