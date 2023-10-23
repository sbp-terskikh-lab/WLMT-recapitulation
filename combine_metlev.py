import pandas as pd
import numpy as np
import os
import dask.dataframe as dd
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10102023"
# results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10122023"

# list_of_samples_of_interest contains paths of the *.metcov.gz files
names = pd.read_csv(
    os.path.join(results_folder, "metcov_samples_of_interest.csv"), header=None
)

# FROM WLMT PAPER -- NOT USED
# all_suitable_sites - a dataframe where column #2 is a list of all sites which pass coverage requirement
selected_sites = pd.read_csv(
    os.path.join(results_folder, "sites_with_good_coverage.csv")
)
selected_sites.set_index("index", inplace=True)

os.makedirs(os.path.join(results_folder, "metlev"), exist_ok=True)


def process_dataframe(name, file):
    data = pd.read_csv(os.path.join(results_folder, "metcov", file))[
        ["index", "metlev"]
    ]
    data.set_index("index", inplace=True)
    data.rename(
        columns={"metlev": f"{name}_metlev"},
        inplace=True,
    )
    data.drop(data.index.values[data.index.str.contains("Y|X|MT")], inplace=True)

    autosome_selected_sites = np.intersect1d(selected_sites.index, data.index)

    metlev_selected_sites = data.loc[autosome_selected_sites, :]

    data.reset_index().to_csv(
        os.path.join(results_folder, "metlev", f"{name}_metlev_full_autosome.csv.gz"),
        compression="gzip",
        index=False,
    )
    metlev_selected_sites.reset_index().to_csv(
        os.path.join(results_folder, "metlev", f"{name}_metlev_selected_sites.csv.gz"),
        compression="gzip",
        index=False,
    )
    return "Finished processing " + name + "\n"
    # return dd.from_pandas(data, chunksize=1000000)


with ProcessPoolExecutor(max_workers=120) as executor:
    futures = [
        executor.submit(process_dataframe, name, file) for name, file in names.values
    ]
    for future in as_completed(futures):
        print(future.result())
    executor.shutdown(wait=True)

# dataframes = []
# with ProcessPoolExecutor(max_workers=120) as executor:
#     futures = [
#         executor.submit(process_dataframe, name, file) for name, file in names.values
#     ]
#     for future in as_completed(futures):
#         dataframes.append(future.result())
#     executor.shutdown(wait=True)
# data = pd.concat(dataframes, join="outer", axis=1, sort=False)
