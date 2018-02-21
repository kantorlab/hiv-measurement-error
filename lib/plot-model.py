import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys
from itertools import product
from sklearn.metrics import roc_curve, roc_auc_score

datasets = ("5VM", "PL11", "PL19")
models = list(product(datasets, repeat=2))

def plot_roc(model):
    y = pd.read_csv("scratch/model.{}.{}.prediction.csv".format(*model))
    fpr, tpr, _ = roc_curve(y["true"], y["pred"])
    auc = roc_auc_score(y["true"], y["pred"])
    sample = np.arange(0, len(fpr), 1000)
    plt.fill_between(fpr[sample], 0, tpr[sample], color="lightgray")
    plt.plot(fpr[sample], tpr[sample])
    plt.annotate(
        "area = %.3f" % auc,
        xy=(0.5, 0.5),
        xytext=(0.5, 0.5),
        va="center",
        ha="center",
        fontsize=18
    )

sns.set_style("darkgrid")
plt.figure(figsize=(18,18))

for i, model in enumerate(models):
    plt.subplot(3, 3, i+1)
    if (i < 3): plt.title(model[1], fontsize=20)
    if (i % 3 == 0): plt.ylabel(model[0], fontsize=20)
    plot_roc(model)

plt.tight_layout()
plt.savefig(sys.argv[1])
