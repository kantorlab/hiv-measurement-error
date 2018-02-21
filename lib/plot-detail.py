import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from Bio import SeqIO
from matplotlib.ticker import FixedLocator

margin = 1.5
freq_threshold = 0.1

dataset = sys.argv[1]
methods = sys.argv[2:]

if dataset == "PID":
  xmax = 319
else:
  xmax = 348

ref = [set() for i in range(xmax)]
for seq in SeqIO.parse("data/{}.ref.fa".format(dataset), "fasta"):
    for i in range(xmax):
        ref[i].add(str(seq.seq[3*i:3*(i+1)]).upper())

codons = {}
coverage = {}
norm_coverage = {}

for method in methods:
    codons[method] = pd.read_csv("results/{}.{}.codons.csv".format(dataset, method),
                                 index_col="AA_Index",
                                 nrows=xmax)
    coverage[method] = codons[method].sum(axis=1)
    norm_coverage[method] = 100.0 / coverage[method]

ymax = max(coverage[method].max() for method in methods)
ylim = (0, ymax)

fig, ax = plt.subplots(nrows=2*len(methods), sharex=True, figsize=(15,6*len(methods)),
                       gridspec_kw={"height_ratios": [1,4]*len(methods)})

ax[0].set_title(dataset, fontsize=24)

xticks = [1] + list(range(10, xmax+1, 10))
xlabels = [1] + list(range(10, 100, 10)) + [1] + list(range(10, xmax-99, 10))
ax[-1].set_xlim(1-margin, xmax+margin)
ax[-1].set_xticks(xticks)
ax[-1].set_xticklabels(xlabels, fontsize=9, rotation=90)
ax[-1].set_xlabel("Codon Position", size=14)
ax[-1].xaxis.set_minor_locator(FixedLocator(list(range(1, xmax+1))))
ax[-1].xaxis.set_tick_params(direction="out", which="both", top="off")

for i, method in enumerate(methods):

    ax1 = ax[2*i]
    ax2 = ax[2*i+1]

    ax1.set_ylim(ylim)
    ax1.set_yticks(ylim)
    ax1.set_yticklabels(list(map("{:,d}".format, ylim)))
    ax1.set_ylabel("Coverage", size=14)

    ax1.bar(np.arange(0.5, xmax+0.5), coverage[method], width=1.0, color="k",
            edgecolor="k")

    ax1.get_xaxis().set_tick_params(direction="out", which="both", top="off")
    ax1.text(1, -0.2*ymax, "Protease", fontsize="14")
    ax1.text(100, -0.2*ymax, "Reverse Transcriptase", fontsize="14")

    axt = ax2.twinx()
    axt.set_ylabel(method, size=32, rotation=270)
    axt.yaxis.set_major_locator(plt.NullLocator())

    ax2.set_ylabel("Codon Frequency", size=14)
    ax2.set_yscale("log")
    ax2.set_ylim(0.09,110)
    yticks = [0.1, 1, 5, 10, 20, 50, 100]
    ax2.set_yticks(yticks)
    ax2.set_yticklabels(['%g%%' % i for i in yticks], size=14)
    ax2.get_yaxis().set_tick_params(direction="out", which="both", top="off")
    ax2.grid(True, which="major")
    ax2.axvline(x=100, color="k")

    variants = ([], [], [])
    for j, row in codons[method].iterrows():
        j = int(j)
        norm = norm_coverage[method][j]
        for codon in row.index:
            freq = norm * row[codon]
            if freq > freq_threshold:
                variants[0].append(j+1)
                variants[1].append(freq)
                if codon in ref[j]:
                    variants[2].append("r")
                else:
                    variants[2].append("k")

    ax2.scatter(variants[0], variants[1], color=variants[2], alpha=0.5,
                marker=".", facecolor="none")

fig.tight_layout()
fig.savefig("results/{}-detail.eps".format(dataset))
