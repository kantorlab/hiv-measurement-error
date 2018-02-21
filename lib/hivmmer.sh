#!/bin/bash
set -e
hivmmer --id $1 --fq1 $1_1.fastq --fq2 $1_2.fastq --ref $2.hmm --cpu $3
mv $1.hmmsearch2.codons.csv $4
