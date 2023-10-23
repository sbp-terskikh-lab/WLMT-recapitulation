import re
import os
import pandas as pd

# results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10122023"
results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10102023"

data_dir = os.path.join(results_folder, "methylation")
files = pd.Series(os.listdir(data_dir))
files = files[files.str.endswith(".bismark.cov.gz")]
names = files.str.split("_").apply(lambda x: x[0])
locations_of_samples_of_interest = pd.concat([names, files], axis=1)
locations_of_samples_of_interest.to_csv(
    f"{results_folder}/cov_samples_of_interest.csv",
    header=False,
    index=False,
    sep=",",
)
