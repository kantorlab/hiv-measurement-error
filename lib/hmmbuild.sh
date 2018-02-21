#!/bin/bash
set -e
FA1=$1
FA2=$2.fa
HMM=$2.hmm
python lib/trim-reference.py $FA1 >$FA2
rm -f $HMM
hmmbuild $HMM $FA2
