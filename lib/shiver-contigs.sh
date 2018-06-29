#!/bin/bash
DIR=$PWD
mkdir -p $1
cd $1
bash $DIR/$2 $DIR/${3%.log} $DIR/$4 $DIR/$5 $6
