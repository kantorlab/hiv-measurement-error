import numpy as np
import pandas as pd
import pysam
import sys
from itertools import product

nts = ("A", "C", "G", "T")
columns = list(map("".join, product(list(map("".join, product(nts, nts))), nts))) + ["del", "ins"]

n = 384
table = pd.DataFrame(data=0, index=np.arange(n), columns=columns, dtype=int)

columns = frozenset(columns)

nreads = 0
npassed = 0
for sam in sys.argv[1:]:
  for aln in pysam.Samfile(sam, "rb").fetch(until_eof=True):
    nreads += 1
    if not (aln.is_unmapped or 'I' in aln.cigar or 'D' in aln.cigar):
      npassed += 1
      frame_start = 3 * int((aln.reference_start + 2) / 3)
      qstart = aln.qstart + (frame_start - aln.reference_start)
      for i in range(qstart, aln.qend, 3):
        codon = aln.seq[i:i+3]
        j = int((frame_start + i - qstart) / 3)
        if j < n and codon in columns:
          table[codon][j] += 1

print(nreads, "reads", file=sys.stderr)
print(npassed, "passed", file=sys.stderr)

table.to_csv(sys.stdout, index_label="AA_Index")
