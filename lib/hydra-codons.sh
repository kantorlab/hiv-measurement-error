#!/bin/bash
set -e

PY=$1
BAM1=$2
BAM2=$3
OUT=$4
TMP=$2.csv

python $PY $BAM1 $BAM2 >$TMP
mv $TMP $OUT
