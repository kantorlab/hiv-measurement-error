import pandas as pd
import numpy as np
import sys
from itertools import combinations
from scipy.stats import friedmanchisquare, wilcoxon
variants = pd.read_csv(sys.argv[1])
variants = variants[variants.reference == 1]
methods = sys.argv[2].split(",")
print(friedmanchisquare(*[variants[method] for method in methods]))
pvalues = []
pairs = []
for m1, m2 in combinations(methods, 2):
    t = wilcoxon(variants[m1], variants[m2])
    print(m1, m2, t)
    pvalues.append(t.pvalue)
    pairs.append("{} {}".format(m1, m2))
print(pvalues)
n = len(pvalues)
pvalues = np.array(pvalues)
order = np.argsort(pvalues)
holm = np.array([np.min([1, i]) for i in (pvalues[order] * (n - np.arange(1, n+1) + 1))])
print(order)
print(np.argsort(order))
for pair, p in zip(pairs, holm[np.argsort(order)]):
    print(pair, p)
