import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

ticks = [
  [0, 0.2, 0.4, 0.6, 0.8, 1.0],
  [0, 0.5, 1.0],
  [0, 0.1, 0.9, 1.0]
]

datasets = sys.argv[1].split(",")
methods  = sys.argv[2].split(",")
csvfiles = sys.argv[3:-1]
epsfile  = sys.argv[-1]

assert len(datasets) == len(csvfiles)

sns.set_style("darkgrid")
plt.figure(figsize=(6*len(datasets),6))

for i, (dataset, csvfile) in enumerate(zip(datasets, csvfiles)):
  variants = pd.read_csv(csvfile, usecols=["reference"]+methods)
  variants = variants[variants.reference == 1]
  del variants["reference"]
  variants = variants.melt(value_vars=methods)
  plt.subplot(1, len(datasets), i+1)
  plt.title(dataset, fontsize=20)
  plt.yticks(ticks[i])
  plt.ylim([-0.05, 1.05])
  sns.swarmplot(x="variable", y="value", data=variants)
  if i==0:
    plt.ylabel("Frequency", fontsize=16)
  else:
    plt.ylabel("")
  if i==1:
    plt.xlabel("Method", fontsize=16)
  else:
    plt.xlabel("")

plt.tight_layout()
plt.savefig(epsfile)
