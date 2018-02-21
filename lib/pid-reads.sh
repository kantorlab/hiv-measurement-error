#!/bin/bash
set -e
cd scratch
cat PrimerID/scripts/Comparing\ estimators/3*/3*_R1.fastq >PID_1.fastq
cat PrimerID/scripts/Comparing\ estimators/3*/3*_R2.fastq >PID_2.fastq
