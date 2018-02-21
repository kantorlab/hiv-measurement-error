#!/bin/bash
set -e
PY=$1
ID=$4
OUT=$5
FASTA=scratch/${ID}.collapsed.fa
N=$(grep -c '^>' $FASTA)
echo $N
python $PY $N >$OUT
