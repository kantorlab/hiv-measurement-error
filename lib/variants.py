import os
import pandas as pd
import numpy as np
import sys
from Bio import SeqIO
from itertools import product

freq_threshold = 0.01

reffile  = sys.argv[1]
methods  = sys.argv[2].split(",")
csvfiles = sys.argv[3:-1]
outfile  = sys.argv[-1]

assert len(methods) == len(csvfiles)

if os.path.basename(reffile).startswith("PID"):
  xmax = 319
else:
  xmax = 348

nt = ("A", "C", "G", "T")
codons = list(map("".join, product(map("".join, product(nt, nt)), nt)))

index = []
for i in range(xmax):
  for codon in codons:
    index.append((i, codon))

variants = pd.DataFrame(0,
                        columns=["reference"] + methods,
                        index=pd.MultiIndex.from_tuples(index, names=["pos", "codon"]))

# Mark correct codons from reference
for seq in SeqIO.parse(reffile, "fasta"):
  for i in range(xmax):
    codon = str(seq.seq[3*i:3*(i+1)]).upper()
    if codon in codons:
        variants.loc[(i, codon), "reference"] = 1

for method, csvfile in zip(methods, csvfiles):
  codons = pd.read_csv(csvfile)
  codons = codons.div(codons.sum(axis=1), axis=0)
  codons[codons < freq_threshold] = 0
  codons = codons.div(codons.sum(axis=1), axis=0)
  for i, codon in variants.index:
    variants.loc[(i, codon), method] = codons.loc[i, codon]

variants.fillna(0).to_csv(outfile)
