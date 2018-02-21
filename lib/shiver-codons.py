import numpy as np
import pandas as pd
import pysam
import sys
from Bio import AlignIO
from itertools import product

nts = ("A", "C", "G", "T")
columns = list(map("".join, product(list(map("".join, product(nts, nts))), nts))) + ["del", "ins"]

msa = AlignIO.read(sys.argv[1], "stockholm")
assert msa[0].name == "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
hxb2 = {}
i = 0
j = 0
for b1, b2 in zip(msa[0].seq, msa[-1].seq):
  if b1 != "-" and b2 != "_":
    hxb2[j] = i
  if b1 != "-":
    i += 1
  if b2 != "-":
    j += 1

n = int((i+2)/3)

table = pd.DataFrame(data=0, index=np.arange(n), columns=columns, dtype=int)

columns = frozenset(columns)

nreads = 0
npassed = 0
for aln in pysam.Samfile(sys.argv[2], "rb").fetch(until_eof=True):
  nreads += 1
  if not (aln.is_unmapped or 'I' in aln.cigar or 'D' in aln.cigar):
    npassed += 1
    reference_start = hxb2.get(aln.reference_start)
    if reference_start is not None:
      frame_start = 3 * int((reference_start + 2) / 3)
      qstart = aln.qstart + (frame_start - reference_start)
      for i in range(qstart, aln.qend, 3):
        codon = aln.seq[i:i+3]
        j = int((frame_start + i - qstart) / 3)
        if j < n and codon in columns:
          table[codon][j] += 1

print(nreads, "reads", file=sys.stderr)
print(npassed, "passed", file=sys.stderr)

table.to_csv(sys.stdout, index_label="AA_Index")
