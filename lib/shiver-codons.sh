#!/bin/bash
set -e

ID=$1
OUT=$2
DIR=scratch/$ID.shiver
REF=$DIR/ref.pol
STO=$DIR/$ID.sto
TMP=$DIR/codons.csv

hmmbuild $REF.hmm $REF.fa
hmmalign --mapali $REF.fa $REF.hmm $DIR/${ID}_remap_ref.fasta >$STO

python lib/shiver-codons.py $STO $DIR/${ID}_remap.bam >$TMP

mv $TMP $OUT
