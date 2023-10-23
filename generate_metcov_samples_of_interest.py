import re
import os
import pandas as pd


results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10102023"
# results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10122023"

data_dir = os.path.join(results_folder, "metcov")
data = os.listdir(data_dir)
files = pd.Series(data)
locations_of_samples_of_interest = files[files.str.endswith(".metcov.gz")]
names = files.str.split(".").apply(lambda x: x[0])
locations_of_samples_of_interest = pd.concat([names, files], axis=1)
locations_of_samples_of_interest.name = None
locations_of_samples_of_interest.to_csv(
    f"{results_folder}/metcov_samples_of_interest.csv",
    header=False,
    index=False,
    sep=",",
)
