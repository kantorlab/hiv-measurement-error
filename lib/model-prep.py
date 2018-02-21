import pandas as pd
import numpy as np
import sys

from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO
from itertools import product

min_score = 1.5

nts = ("A", "C", "G", "T")
codons = list(map("".join, product(list(map("".join, product(nts, nts))), nts))) + ["del", "ins"]
lookup = dict((codon, i) for (i, codon) in enumerate(codons))

def codon_table(reads, nreads, refcodons, hmmer):

  fields = [
    "error", "codon",
    "readpos", "readlen", "readfreq",
    "alnpos", "alnlen", "alnprob", "alnscore",
    "revcomp", "frame", "pos"]
  print(",".join(fields))

  for hit in hmmer.hits:

    row = {}

    id, _, frame = hit.id.rpartition("-")
    count = int(id.partition("-")[2])
    row["readfreq"] = count*nreads

    if frame.endswith("'"):
      row["revcomp"] = 1
      seq = reads[id].seq.reverse_complement()
      offset = int(frame[:-1])
    else:
      row["revcomp"] = 0
      seq = reads[id].seq
      offset = int(frame)
    row["readlen"] = len(seq)
    row["frame"] = offset

    for hsp in hit.hsps:

      if (hsp.bitscore / hsp.hit_span) < min_score: continue

      row["alnscore"] = hsp.bitscore
      row["alnlen"] = len(hsp.aln[1].seq)

      i = 3*hsp.hit_start + offset
      j = 0

      for aa, pp in zip(hsp.aln[1].seq, hsp.aln_annotation["PP"]):
        if aa.islower():
          codon = "ins"
        elif aa == "-":
          codon = "del"
        else:
          codon = str(seq[i:i+3])
          assert aa == Seq.translate(codon)
        row["error"] = int(codon not in refcodons[hsp.query_start+j])
        row["codon"] = codon
        row["readpos"] = float(i) / row["readlen"]
        row["alnpos"] = float(j) / row["alnlen"]
        row["alnprob"] = pp.replace("*", "10").replace(".", "-1")
        row["pos"] = hsp.query_start+j
        print(",".join(map(str, list(map(row.get, fields)))))
        j += (codon != "ins")
        i += 3*(codon != "del")


if __name__ == "__main__":

  reads = SeqIO.index(sys.argv[1], "fasta")
  print("indexed reads", file=sys.stderr)

  nreads = 1.0 / float(sys.argv[2])

  refs = SeqIO.index(sys.argv[3], "fasta")
  print("indexed references", file=sys.stderr)

  n = max(len(ref.seq) for ref in list(refs.values()))
  refcodons = [set() for x in range(n)]
  for ref in refs.values():
    for i in range(0, len(ref.seq), 3):
      j = int(i/3)
      refcodons[j].add(str(ref.seq[i:i+3]).upper())

  hmmer = SearchIO.read(sys.argv[4], "hmmer3-text")
  print("parsed hmmer", file=sys.stderr)

  codon_table(reads, nreads, refcodons, hmmer)

# vim: expandtab sw=2 ts=2
