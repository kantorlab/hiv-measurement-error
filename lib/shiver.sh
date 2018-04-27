#!/bin/bash
set -e

ID=$1
OUT=scratch/shiver.${ID}
TMP=scratch/shiver.${ID}.tmp

rm -rf $TMP $OUT
mkdir -p $TMP
cd $TMP
ln -s ../shiver/tools .
ln -s ../shiver/info .
ln -s ../shiver/*.sh .
ln -s ../shiver/*.fa .

./shiver_init.sh reference config.sh ref.pol.fa adapters.fa primers.fa

./shiver_align_contigs.sh reference config.sh ../iva.${ID}/contigs.fasta $ID

./shiver_map_reads.sh reference config.sh ../iva.${ID}/contigs.fasta $ID \
  ${ID}.blast ${ID}_cut_wRefs.fasta ../${ID}_?.fastq

cd ../../
mv $TMP $OUT
