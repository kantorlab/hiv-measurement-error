#!/bin/bash
set -e

ID=$1
TARGET=scratch/$ID.shiver
TMP=scratch/$ID.shiver.tmp

rm -rf $TMP $TARGET
mkdir -p $TMP
cd $TMP
ln -s ../shiver/tools .
ln -s ../shiver/info .
ln -s ../shiver/*.sh .
ln -s ../shiver/*.fa .

./shiver_init.sh reference config.sh ref.pol.fa adapters.fa primers.fa

./shiver_align_contigs.sh reference config.sh ../${ID}.iva/contigs.fasta $ID

./shiver_map_reads.sh reference config.sh ../${ID}.iva/contigs.fasta $ID \
  ${ID}.blast ${ID}_cut_wRefs.fasta ../${ID}_?.fastq

cd ../../
mv $TMP $TARGET
