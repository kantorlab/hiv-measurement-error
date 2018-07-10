#!/usr/bin/env python
import sys
from Bio import SearchIO

min_score = 1.5

if __name__ == "__main__":

  hmmer = SearchIO.read(sys.argv[1], "hmmer3-text")
  passed = {}

  for hit in hmmer.hits:

    id, _, frame = hit.id.rpartition("-")
    count = int(id.partition("-")[2])

    for hsp in hit.hsps:

      if (hsp.bitscore / hsp.hit_span) >= min_score:
        passed[id] = count

  print(sum(passed.values()))

# vim: expandtab sw=2 ts=2
