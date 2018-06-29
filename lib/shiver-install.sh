#!/bin/bash
set -e
ZIP=$1
DIR=${ZIP%.zip}
rm -rf $DIR scratch/shiver
unzip $ZIP -d scratch
mv $DIR scratch/shiver
sed -i~ -e 's|#!/usr/bin/env python2|#!/usr/bin/env python|' scratch/shiver/tools/*.py
2to3 -w scratch/shiver/tools
