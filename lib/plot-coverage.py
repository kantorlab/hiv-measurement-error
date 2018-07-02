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

re_cigar = re.compile(r'(\d+)([A-Z]{1})')
fragment_types = ["Replicate", "Overlapping", "Paired"]
hxb2_len = 9719
xmax = 1200
ops = {"M": 1, "D": 1, "I": 0}

def cigar2len(cigar):
  l = 0
  for n, op in re_cigar.findall(cigar):
      l += ops[op] * int(n)
  return l

raw = pd.read_csv(tsvfile, delimiter="\t", names=["start", "cigar", "fragment"])

raw["length"] = raw["cigar"].apply(cigar2len)
raw["fragment"] = np.abs(raw["fragment"])

fragments = np.zeros((2**14, 3))
coverage = np.zeros((2**14, 3))
for row in raw.itertuples():
  i = min(2, int(row.fragment/250))
  fragments[int(row.fragment/10),i] += 1
  coverage[int(row.start/50):int((row.start+row.length)/50),i] += 1
fragments = pd.DataFrame(fragments[:int(xmax/10)+1,:], columns=fragment_types)
coverage = pd.DataFrame(coverage[:int(hxb2_len/50),:], columns=fragment_types)

sns.set_style("darkgrid")

f, (ax0, ax1) = plt.subplots(1, 2, figsize=(12,4), gridspec_kw={"width_ratios":[1, 2]})

fragments.plot(kind="bar", ax=ax0, stacked=True)
ax0.set_title("a")
ax0.set_xlabel("Fragment size")
ax0.set_xticks([0,xmax/10])
ax0.set_xticklabels([0,xmax])
ax0.set_ylabel("Frequency")

coverage.plot(kind="bar", ax=ax1, stacked=True)
ax1.set_title("b")
ax1.set_xlabel("Genome position in HXB2 coordinates")
ax1.set_xticks([0,int(hxb2_len/50)])
ax1.set_xticklabels([0,hxb2_len])
ax1.legend()

f.tight_layout()
f.savefig(epsfile)
