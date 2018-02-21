#!/bin/bash
set -e
cd scratch
rm -f v1.3.3.zip
wget https://github.com/ChrisHIV/shiver/archive/v1.3.3.zip
rm -rf shiver-1.3.3
unzip v1.3.3.zip
rm -f shiver
ln -s shiver-1.3.3 shiver
sed -i~ -e 's|#!/usr/bin/env python2|#!/usr/bin/env python|' shiver/tools/*.py
2to3 -w shiver/tools
cd ..
cp data/*.fa scratch/shiver/
