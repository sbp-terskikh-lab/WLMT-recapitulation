import pandas as pd
import subprocess
import os
from tqdm import tqdm
import dask.dataframe as dd
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10102023"
# results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10122023"

names = pd.read_csv(
    os.path.join(results_folder, "metcov_samples_of_interest.csv"), header=None
)


def process_dataframe(name, file):
    data = pd.read_csv(os.path.join(results_folder, "metcov", file))[
        ["index", "coverage"]
    ]
    data.set_index("index", inplace=True)
    data.rename(
        columns={"coverage": f"{name}_coverage"},
        inplace=True,
    )
    return data
    # return dd.from_pandas(data, chunksize=1000000)


dataframes = []
with ProcessPoolExecutor() as executor:
    futures = [
        executor.submit(process_dataframe, name, file) for name, file in names.values
    ]
    for future in as_completed(futures):
        dataframes.append(future.result())
    executor.shutdown(wait=True)

data = pd.concat(dataframes, join="inner", axis=1)

d0_covered_5_90 = data[data.ge(5).sum(axis=1) / data.shape[1] >= 0.90]

d0_covered_5_90.reset_index(inplace=True)

d0_covered_5_90.to_csv(
    os.path.join(results_folder, "sites_with_good_coverage.csv"),
    index=False,
)
