#!/bin/bash
set -e

ID=$1
OUT=$2
DIR=scratch/shiver.$ID
REF=$DIR/ref.pol
STO=$DIR/${ID}.sto
BAM=$DIR/${ID}_remap.bam

hmmbuild $REF.hmm $REF.fa
hmmalign --mapali $REF.fa $REF.hmm $DIR/${ID}_remap_ref.fasta >$STO

python lib/pileup.py $STO $BAM
mv ${BAM}.csv $OUT
