import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

from Bio import SeqIO
from itertools import groupby

datasets = sys.argv[1].split(",")
methods  = sys.argv[2].split(",")
reffiles = sys.argv[3:3+len(datasets)]
csvfiles = sys.argv[3+len(datasets):-1]
epsfile  = sys.argv[-1]

assert len(csvfiles) == len(datasets)*len(methods)

freq_threshold = 0.001
xmax = 348

def load_data(reffile, csvfiles):
  ref = [set() for i in range(xmax)]
  for seq in SeqIO.parse(reffile, "fasta"):
    for i in range(xmax):
      ref[i].add(str(seq.seq[3*i:3*(i+1)]).upper())

  x = {}
  y = {}

  for method, f in zip(methods, csvfiles):
    codons = pd.read_csv(f, nrows=xmax)
    errors = []
    for j, row in codons.iterrows():
      j = int(j)
      row = row / row.sum()
      for codon in row.index:
        if row[codon] >= freq_threshold and (not codon in ref[j]):
          errors.append(row[codon])

    errors = sorted(-1*np.log10(errors))
    count = 0
    x[method] = []
    y[method] = []
    for e, group in groupby(errors):
      count += len(list(group))
      x[method].append(e)
      y[method].append(count)

  return x, y

sns.set_style("darkgrid")
plt.figure(figsize=(6*len(datasets),6))

for i, dataset in enumerate(datasets):
  plt.subplot(1, len(datasets), i+1)
  plt.title(dataset, fontsize=20)
  plt.xlim([-np.log10(0.02), 3])
  plt.xticks(
      [1, -np.log10(0.05), -np.log10(0.02), 2, -np.log10(0.005), -np.log10(0.0025), 3],
      ["10%", "5%", "2%", "1%", "0.5%", "0.25%", "0.1%"])
  plt.ylim((-50, 1600))
  plt.axvline(2, c="k", lw=0.5)
  x, y = load_data(reffiles[i], csvfiles[i::len(datasets)])
  for method in methods:
    plt.plot(x[method], y[method], label=method)
    plt.scatter(x[method][0], y[method][0], label=None)
  if i==0:
    plt.ylabel("# of Erroneous Variants", fontsize=16)
    plt.legend(loc="upper left")
  if i==1:
    plt.xlabel("Frequency Threshold (log-scaled)", fontsize=16)

plt.tight_layout()
plt.savefig(epsfile)
