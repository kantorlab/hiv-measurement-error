#!/bin/bash
set -e

ID=$1
REF=$2
CPU=$3
IDX=scratch/${ID}.coverage.idx

rm -f $IDX.*
bowtie2-build $REF $IDX
bowtie2 --very-sensitive -X 10000 --threads $CPU -x $IDX -1 scratch/${ID}_1.fastq -2 scratch/${ID}_2.fastq | samtools view -b -S - >scratch/${ID}.coverage.bam
samtools view -f 3 scratch/${ID}.coverage.bam | cut -f 4,6,9 >scratch/${ID}.coverage.tsv
