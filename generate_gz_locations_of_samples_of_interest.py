import re
import os
import pandas as pd

files = pd.Series(os.listdir("./data/GSE121141_RAW"))
files = files[files.str.endswith(".gz")]
names = files.str.split(".").apply(lambda x: x[0])
locations_of_samples_of_interest = pd.concat([names, files], axis=1)
locations_of_samples_of_interest.to_csv(
    "./data/locations_of_samples_of_interest.csv", index=False, sep=','
)
