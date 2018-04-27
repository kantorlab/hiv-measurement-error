#!/bin/bash
set -e

PY=$1
DIR=$2
OUT=$3

python $PY $DIR/align.bam >${OUT}~
mv ${OUT}~ $OUT
