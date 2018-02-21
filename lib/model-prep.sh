#!/bin/bash
set -e
PY=$1
HMM=$2
ID=$3
REF=data/${ID}.ref.fa
FASTA=scratch/${ID}.collapsed.fa
N=$(grep -c '^>' $FASTA)
echo $N
python $PY $FASTA $N $REF scratch/${ID}.${HMM}.txt >scratch/model.${ID}.csv
