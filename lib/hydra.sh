#!/bin/bash
set -e
ID=$1
TMP=scratch/hydra.${ID}.tmp
OUT=scratch/hydra.${ID}
rm -rf $TMP $OUT
quasitools hydra -o $TMP scratch/${ID}_1.fastq scratch/${ID}_2.fastq
mv $TMP $OUT
