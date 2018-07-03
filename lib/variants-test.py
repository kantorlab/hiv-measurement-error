import pandas as pd
import sys
from itertools import combinations
from scipy.stats import friedmanchisquare, wilcoxon
variants = pd.read_csv(sys.argv[1])
variants = variants[variants.reference == 1]
methods = sys.argv[2].split(",")
print(friedmanchisquare(*[variants[method] for method in methods]))
for m1, m2 in combinations(methods, 2):
    print(m1, m2, wilcoxon(variants[m1], variants[m2]))
