import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

from Bio import SeqIO
from itertools import groupby

freq_threshold = 0.001
xmax = 319

reffile  = sys.argv[1]
varfile  = sys.argv[2]
csvfiles = sys.argv[3:-1]
epsfile  = sys.argv[-1]

methods = [f.split(".")[1] for f in csvfiles]

def load_data():
  seqs = list(SeqIO.parse(reffile, "fasta"))
  reference = [set() for i in range(xmax)]
  for seq in seqs:
    for i in range(xmax):
      reference[i].add(str(seq.seq[3*i:3*(i+1)]).upper())

  x = {}
  y = {}
  for method, csvfile in zip(methods, csvfiles):
    codons = pd.read_csv(csvfile, nrows=xmax)
    errors = []
    for j, row in codons.iterrows():
      if row.sum() > 0 and j != 236:
        row = row / row.sum()
        for codon in row.index:
          if row[codon] >= freq_threshold and (not codon in reference[j]):
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
plt.figure(figsize=(12,6))

plt.subplot(1, 2, 1)
plt.title("Cumulative Errors", fontsize=20)

plt.xlim([-np.log10(0.05), 3])
plt.xticks(
  [-np.log10(0.05), -np.log10(0.02), 2, -np.log10(0.005), -np.log10(0.0025), 3],
  ["5%", "2%", "1%", "0.5%", "0.25%", "0.1%"])
plt.ylim((-25, 700))
x, y = load_data()
for method in methods:
  plt.plot(x[method], y[method], label=method)
  plt.scatter(x[method][0], y[method][0], label=None)
plt.ylabel("# of Erroneous Variants", fontsize=16)
plt.legend(loc="upper left")
plt.xlabel("Frequency Threshold (log-scaled)", fontsize=16)

plt.subplot(1, 2, 2)
plt.title("Variant Distributions", fontsize=20)

variants = pd.read_csv(varfile, usecols=["reference"]+methods)
variants = variants[variants.reference == 1]
del variants["reference"]
variants = variants.melt(value_vars=methods)
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.ylim([-0.05, 1.05])
sns.swarmplot(x="variable", y="value", data=variants)
plt.ylabel("Frequency", fontsize=16)
plt.xlabel("Method", fontsize=16)

plt.tight_layout()
plt.savefig(epsfile)
