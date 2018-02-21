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

methods = ["hydra", "hivmmer", "hivmmer-ml", "pidalyse"]

def load_data():
  seqs = list(SeqIO.parse("data/PID.ref.fa", "fasta"))
  reference = [set() for i in range(xmax)]
  for seq in seqs:
    for i in range(xmax):
      reference[i].add(str(seq.seq[3*i:3*(i+1)]).upper())

  x = {}
  y = {}
  for method in methods:
    codons = pd.read_csv("results/PID.{}.codons.csv".format(method),
                         index_col="AA_Index",
                         nrows=xmax)
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

plt.xlim([-np.log10(0.02), 3])
plt.xticks(
  [-np.log10(0.02), 2, -np.log10(0.005), -np.log10(0.0025), 3],
  ["2%", "1%", "0.5%", "0.25%", "0.1%"])
plt.ylim((-25, 500))
x, y = load_data()
for method in methods:
  plt.plot(x[method], y[method], label=method)
  plt.scatter(x[method][0], y[method][0], label=None)
plt.ylabel("# of Erroneous Calls", fontsize=16)
plt.legend(loc="upper left")
plt.xlabel("Frequency Threshold (log-scaled)", fontsize=16)

plt.subplot(1, 2, 2)
plt.title("Calls", fontsize=20)

calls = pd.read_csv("results/PID.calls.csv",
                    usecols=methods,
                    nrows=xmax)
calls = calls.melt(value_vars=methods)
plt.xlabel("")
plt.ylabel("")
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.ylim([-0.05, 1.05])
sns.swarmplot(x="variable", y="value", data=calls)
plt.ylabel("Frequency", fontsize=16)
plt.xlabel("Method", fontsize=16)

plt.tight_layout()
plt.savefig(sys.argv[1])
