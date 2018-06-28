#!/bin/bash
set -e
cd scratch
rm -rf shiver shiver-1.3.3
unzip shiver-1.3.3.zip
ln -s shiver-1.3.3 shiver
sed -i~ -e 's|#!/usr/bin/env python2|#!/usr/bin/env python|' shiver/tools/*.py
2to3 -w shiver/tools
cd ..
cp data/*.fa scratch/shiver/
