import csv
import pysam
import sys
from Bio import SeqIO

ref = str(next(SeqIO.parse(sys.argv[1], "fasta")).seq).upper()
n = len(ref)

pol_start = 2253

nreads = 0
npassed = 0

writer = csv.writer(sys.stdout)
writer.writerow(("error", "qual", "readpos"))

for aln in pysam.Samfile(sys.argv[2], "rb").fetch(until_eof=True):
  nreads += 1
  if not aln.is_unmapped:
    npassed += 1
    for i, (nt, qual) in enumerate(zip(aln.seq, aln.qual)):
      j = aln.reference_start + i - pol_start + 1
      if j > 0 and j < n and ref[j] != "-":
          writer.writerow((int(nt != ref[j]), ord(qual)-33, i))
    if nreads > 10: break

# vim: expandtab sw=2 ts=2
