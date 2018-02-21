import numpy as np
import pandas as pd
import pickle
import sys
from itertools import product
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

cpu = int(sys.argv[1])
np.random.seed(10975791)

datasets = ("5VM", "PL11", "PL19")
data = {}

# load data
for dataset in datasets:
    data[dataset] = pd.read_csv(
            "scratch/model.{}.csv".format(dataset),
            converters={"alnprob": lambda x: np.uint8(x.replace(".", "-1").replace("*", "10"))}
            )
    del data[dataset]["codon"]
    del data[dataset]["pos"]

def save_prediction(y_true, y_pred, outname):
    df = pd.DataFrame(0, index=np.arange(len(y_true)), columns=["true", "pred"])
    df["true"] = y_true
    df["pred"] = y_pred
    df.to_csv(outname, index=False)

def save_importance(names, values, outname):
    df = pd.DataFrame(0, index=np.arange(len(names)), columns=["Feature", "Importance"])
    df["Feature"] = names
    df["Importance"] = values
    df.to_csv(outname, index=False)

# run a model for each pair of training/validation data
models = list(product(datasets, repeat=2))

for training, validation in models:
    model = XGBClassifier(max_depth=8, n_jobs=cpu, random_state=np.random.randint(2**31))

    if training == validation:
        y_train, y_validate, X_train, X_validate = train_test_split(
            data[training].iloc[:,0].as_matrix(),
            data[training].iloc[:,1:].as_matrix(),
            test_size=0.5,
            random_state=857719029
        )
        print(y_train)
        print(y_validate)
    else:
        y_train = data[training].iloc[:,0]
        X_train = data[training].iloc[:,1:]
        y_validate = data[validation].iloc[:,0]
        X_validate = data[validation].iloc[:,1:]

    model.fit(X_train, y_train)
    pickle.dump(model, open("scratch/model.{}.{}.bin".format(training, validation), "wb"))

    prediction = model.predict_proba(X_validate)[:,1]

    print(confusion_matrix(y_validate, prediction > 0.5))
    save_prediction(
        y_validate,
        prediction,
        "scratch/model.{}.{}.prediction.csv".format(training, validation))
    save_importance(
        data[training].columns[1:],
        model.feature_importances_,
        "scratch/model.{}.{}.importance.csv".format(training, validation))
