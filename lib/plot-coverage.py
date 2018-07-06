import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import seaborn as sns
import sys

tsvfile = sys.argv[1]
epsfile = sys.argv[2]

fragment_types = ["Replicate", "Overlapping", "Paired"]

raw = pd.read_csv(tsvfile, delimiter="\t", names=["start", "fragment"])
raw["fragment"] = np.abs(raw["fragment"])

fragments = np.zeros((2**14, 3))
for row in raw.itertuples():
  i = min(2, int(row.fragment/250))
  fragments[int(row.fragment/10),i] += 1
fragments = pd.DataFrame(fragments[:101,:], columns=fragment_types)

sns.set_style("darkgrid")

fig = plt.figure(figsize=(6, 6))
ax = fig.gca()
fragments.plot(kind="bar", ax=ax, stacked=True)
ax.set_xlabel("Fragment size in bowtie2 alignment", fontsize=14)
ax.set_xticks([0, 25, 50, 100])
ax.set_xticklabels([0, 250, 500, 1000])
ax.set_ylabel("Frequency", fontsize=14)

fig.tight_layout()
fig.savefig(epsfile)
