import os
import pandas as pd
import numpy as np
import sys
from Bio import SeqIO

freq_threshold = 0.01
dataset = sys.argv[1]
methods = sys.argv[2:]

if dataset == "PID":
  xmax = 319
else:
  xmax = 348

reference = [set() for i in range(xmax)]
for seq in SeqIO.parse("data/{}.ref.fa".format(dataset), "fasta"):
  for i in range(xmax):
    reference[i].add(str(seq.seq[3*i:3*(i+1)]).upper())

index = []
for i, codons in enumerate(reference):
  for codon in codons:
    index.append((i, codon))

calls = pd.DataFrame(0,
                     columns=methods,
                     index=pd.MultiIndex.from_tuples(index, names=["pos", "codon"]))

for method in methods:
  codons = pd.read_csv("results/{}.{}.codons.csv".format(dataset, method),
                       index_col="AA_Index")
  if dataset == "PID" and method.startswith("hivmmer"):
    codons.index += 161
  codons = codons.div(codons.sum(axis=1), axis=0)
  codons[codons < freq_threshold] = 0
  codons = codons.div(codons.sum(axis=1), axis=0)
  for i, codon in calls.index:
    calls.loc[(i, codon), method] = codons.loc[i, codon]

calls.fillna(0).to_csv("results/{}.calls.csv".format(dataset))
