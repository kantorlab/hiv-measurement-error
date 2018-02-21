import pandas as pd
import sys
from scipy.stats import kruskal

methods = ["hydra", "shiver", "hivmmer"]

datasets = sys.argv[1:]

for dataset in datasets:
    calls = pd.read_csv(dataset)
    print(dataset, kruskal(*[calls[method] for method in methods]))
