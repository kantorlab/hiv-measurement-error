#!/bin/bash
set -e

cd scratch/PrimerID/scripts/RawData

rm -f *.fastq

ln -s ../../../SRR2097103_1.fastq 3223aR1.fastq
ln -s ../../../SRR2097103_2.fastq 3223aR2.fastq

ln -s ../../../SRR2097104_1.fastq 3223bR1.fastq
ln -s ../../../SRR2097104_2.fastq 3223bR2.fastq

ln -s ../../../SRR2097105_1.fastq 3223cR1.fastq
ln -s ../../../SRR2097105_2.fastq 3223cR2.fastq

ln -s ../../../SRR2097106_1.fastq 3236aR1.fastq
ln -s ../../../SRR2097106_2.fastq 3236aR2.fastq

ln -s ../../../SRR2097107_1.fastq 3236bR1.fastq
ln -s ../../../SRR2097107_2.fastq 3236bR2.fastq

ln -s ../../../SRR2097108_1.fastq 3236cR1.fastq
ln -s ../../../SRR2097108_2.fastq 3236cR2.fastq

