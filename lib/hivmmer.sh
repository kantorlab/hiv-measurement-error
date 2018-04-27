#!/bin/bash
set -e
ID=scratch/$1
REF=$2.hmm
CPU=$3
OUT=$4
hivmmer --id $ID --fq1 ${ID}_1.fastq --fq2 ${ID}_2.fastq --ref $REF --cpu $CPU
mv ${ID}.hmmsearch2.codons.csv $OUT
