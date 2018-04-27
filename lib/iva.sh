#!/bin/bash
set -e
ID=$1
CPU=$2
OUT=scratch/iva.${ID}
TMP=scratch/iva.${ID}.tmp
rm -rf $TMP $OUT
iva -t $CPU -f scratch/${ID}_1.fastq -r scratch/${ID}_2.fastq $TMP
mv $TMP $OUT
