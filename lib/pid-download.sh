#!/bin/bash
set -e

cd scratch

prefetch SRR2097103
prefetch SRR2097104
prefetch SRR2097105
prefetch SRR2097106
prefetch SRR2097107
prefetch SRR2097108

for i in {3..8}
do
	fastq-dump --defline-qual '+' --split-files SRR209710${i}
done

