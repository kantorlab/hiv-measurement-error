#!/bin/bash
set -e
PY=$1
ID=$2
OUT=$3
python $PY scratch/PrimerID/scripts/Alignments/3*/*_QC_1_cons.fasta >scratch/PrimerID/codons.csv
mv scratch/PrimerID/codons.csv $OUT
