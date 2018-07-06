import pandas as pd
import sys

varfile = sys.argv[1]
methods = sys.argv[2].split(",")
outfile = sys.argv[3]

variants = pd.read_csv(varfile)
variants = variants[variants.reference == 0]

top10 = [variants[["pos", "codon", method]].sort_values(method, ascending=False).reset_index(drop=True).head(10)
         for method in methods]

pd.concat(top10, axis=1).to_csv(outfile)
