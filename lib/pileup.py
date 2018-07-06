import numpy as np
import pysam
import sys
from itertools import product

nt = ("A", "C", "G", "T")
codons = list(map("".join, product(map("".join, product(nt, nt)), nt)))
codon_index = dict((codon, i) for (i, codon) in enumerate(codons))

offset = int(sys.argv[1])
n = 348

table = np.zeros((n, len(codons)), dtype=np.uint32)

nreads = 0
npassed = 0
for aln in pysam.Samfile(sys.argv[2], "rb").fetch(until_eof=True):
  nreads += 1
  reference_start = aln.reference_start - offset
  if not (aln.is_unmapped or 'I' in aln.cigar or 'D' in aln.cigar):
    npassed += 1
    reference_start = aln.reference_start - offset
    frame_start = 3 * int((reference_start + 2) / 3)
    qstart = aln.qstart + (frame_start - reference_start)
    for i in range(qstart, aln.qend, 3):
        codon = aln.seq[i:i+3]
        j = int((frame_start + i - qstart) / 3)
        if 0 <= j and j < n and codon in codon_index:
          table[j, codon_index[codon]] += 1

print(nreads, "reads", file=sys.stderr)
print(npassed, "passed", file=sys.stderr)

np.savetxt(sys.argv[3], table, fmt="%g", comments="",
           delimiter=",", header=",".join(codons))
