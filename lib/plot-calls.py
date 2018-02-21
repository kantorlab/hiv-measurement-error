import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys

ticks = [
  [0, 0.2, 0.4, 0.6, 0.8, 1.0],
  [0, 0.1, 0.9, 1.0],
  [0, 0.5, 1.0]
]
methods = ["hydra", "shiver", "hivmmer"]
datasets = ["5VM", "PL1:1", "PL1:9"]

xmax = 348

sns.set_style("darkgrid")
plt.figure(figsize=(6*len(datasets),6))

for i, dataset in enumerate(datasets):
  calls = pd.read_csv("results/{}.calls.csv".format(dataset.replace(":", "")),
                      usecols=methods,
                      nrows=xmax)
  calls = calls.melt(value_vars=methods)
  plt.subplot(1, len(datasets), i+1)
  plt.title(dataset, fontsize=20)
  plt.xlabel("")
  plt.ylabel("")
  plt.yticks(ticks[i])
  plt.ylim([-0.05, 1.05])
  sns.swarmplot(x="variable", y="value", data=calls)
  if i==0:
    plt.ylabel("Frequency", fontsize=16)
  if i==1:
    plt.xlabel("Method", fontsize=16)

plt.tight_layout()
plt.savefig(sys.argv[1])
