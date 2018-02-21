#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys

from Bio import SeqIO
from itertools import product

offset = 463

nts = ("A", "C", "G", "T")
codons = list(map("".join, product(list(map("".join, product(nts, nts))), nts)))
lookup = dict((codon, i) for (i, codon) in enumerate(codons))

def codon_table(fasta):

  nref = 384

  table = np.zeros((nref, len(codons)), dtype=int)

  for f in fasta:
    for seq in SeqIO.parse(f, "fasta"):
      seq = str(seq.seq)
      for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if not ('N' in codon or '-' in codon):
          table[int((offset+i)/3), lookup[codon]] += 1

  return pd.DataFrame(table, index=np.arange(nref), columns=codons)

if __name__ == "__main__":
  codon_table(sys.argv[1:]).to_csv(sys.stdout, index_label="AA_Index")
