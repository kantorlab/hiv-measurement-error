import numpy as np
import pandas as pd
import pickle
import sys
from itertools import product
from xgboost import XGBClassifier

nreads = int(sys.argv[1])

nts = ("A", "C", "G", "T")
columns = list(map("".join, product(list(map("".join, product(nts, nts))), nts))) + ["del", "ins"]

n = 319
table = pd.DataFrame(data=0, index=np.arange(n), columns=columns, dtype=int)

model = pickle.load(open("scratch/model.5VM.PL11.bin", "rb"))

X = pd.read_csv(
        "scratch/model.PID.csv",
        converters={"alnprob": lambda x: np.uint8(x.replace(".", "-1").replace("*", "10"))},
        )

pos = list(X["pos"])
codons = list(X["codon"])
del X["error"]
del X["codon"]
del X["pos"]

prediction = model.predict_proba(X)[:,1]

for j, codon, p, freq in zip(pos, codons, prediction, X["readfreq"]):
    if p < 0.01:
        table[codon][j] += freq * nreads

table.to_csv(sys.stdout, index_label="AA_Index")
