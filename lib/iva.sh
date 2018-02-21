#!/bin/bash
set -e
ID=$1
CPU=$2
TARGET=scratch/$ID.iva
TMP=scratch/$ID.iva.tmp
rm -rf $TARGET $TMP
iva -t $CPU -f scratch/${ID}_1.fastq -r scratch/${ID}_2.fastq $TMP
mv $TMP $TARGET
