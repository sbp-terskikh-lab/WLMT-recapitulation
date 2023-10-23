import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from tqdm import tqdm
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import pickle
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# from sklearn.externals import joblib
import joblib

from scipy.sparse import csc_matrix, coo_matrix, hstack

from scipy.stats import pearsonr, spearmanr

import tables

import time

# Initialize WD
root_folder = "/cjc/miel/martin/WLMT-recapitulation"
results_folder = f"{root_folder}/data/results/results_10102023"
# results_folder = f"{root_folder}/data/results/results_10122023"
output_folder = f"{results_folder}/analysis_result"
os.makedirs(output_folder, exist_ok=True)

pca2 = joblib.load(f"{output_folder}/pca_0.01propnan_filtered.joblib")
prop_nan2 = pd.read_csv(f"{output_folder}/prop_nan.csv")
# np.save(
#     f"{output_folder}/pca_explained_variance_ratio_propnan_filtered.npy",
#     pca2.explained_variance_ratio_,
# )

metlev_merge = pd.read_hdf(f"{output_folder}/metlev_merge.hdf5", key="df")


sample_genes = joblib.load(f"{output_folder}/sample_genes.joblib")


def merge_arrays(arrays):
    """
    Merge multiple NumPy arrays by taking the intersection of their values.
    """
    # Convert arrays to sets
    arrays = [set(array) for array in arrays]
    # Find the intersection of the arrays
    intersection = arrays[0]
    intersection.intersection_update(*arrays[1:])
    # Convert the intersection to a NumPy array
    return np.array(list(intersection))


def merge_arrays_parallel(arrays, num_workers=4):
    """
    Merge multiple NumPy arrays in parallel using a ProcessPoolExecutor.
    """
    start = time.time()
    runtime = 0
    timeout = 60 * 60 * 1  # 1 hours in seconds
    while len(arrays) > 2 and runtime < timeout:
        # Split the arrays into chunks for each worker
        chunk_size = np.ceil(len(arrays) / num_workers).astype(int)
        chunks = [arrays[i : i + chunk_size] for i in range(0, len(arrays), chunk_size)]
        # Merge the arrays in parallel using a ProcessPoolExecutor
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(merge_arrays, chunk) for chunk in chunks]
            arrays = [future.result() for future in futures]
        # Update the runtime
        runtime = time.time() - start
    if runtime >= timeout:
        print("Warning: merge_arrays_parallel timed out")
        return None
    # Merge the results from each worker
    return merge_arrays(arrays)


# merged_sample_genes2 = merge_arrays_parallel(list(sample_genes.values()), num_workers=4)
# np.save(f"{output_folder}/merged_sample_genes.npy", merged_sample_genes)
# propnan_thresh = 0.05
# prop_nan = (metlev_merge == 0).sum(axis=1) / metlev_merge.shape[1]
# threshres = prop_nan < propnan_thresh
# threshres.sum()
# all_meta = pd.read_csv(f"{output_folder}/all_meta.csv")


# a = meer_clock["Weight"].values
# l2a = np.linalg.norm(a)
# mmetage = metlev_merge.loc[meer_clock["index"], :].apply(
#     lambda c: np.dot(c, a) / l2a, axis=0
# )


propnan_thresh = 0.01
suffix = f"{propnan_thresh}propnan_filtered"
prop_nan = (metlev_merge == 0).sum(axis=1) / metlev_merge.shape[1]
pca_data2 = pca_data.loc[:, prop_nan <= propnan_thresh]
pca_data2.loc[:, meer_clock["index"]] = pca_data.loc[:, meer_clock["index"]]
prop_nan.rename("prop_nan")
prop_nan.reset_index().to_csv(f"{output_folder}/prop_nan_{suffix}.csv", index=False)
# Perform PCA on metcov_data
pca2 = PCA(n_components=pca_data2.shape[0])
pca_result2 = pca2.fit_transform(pca_data2)
pca_result2 = pd.DataFrame(
    pca_result2,
    index=pca_data2.index,
    columns=["PC" + str(i) for i in range(1, pca_result2.shape[1] + 1)],
)
joblib.dump(pca2, f"{output_folder}/pca_{suffix}.joblib")
np.save(
    f"{output_folder}/pca_explained_variance_ratio_{suffix}.npy",
    pca2.explained_variance_ratio_,
)

# with open(f"{output_folder}/pca_{propnan_thresh}propnan_filtered.pickle", "wb") as f:
#     pickle.dump(pca2, f)


pca_wmeta2 = pd.merge(
    metage_wmeta, pca_result2, left_on="run", right_index=True, how="inner"
)
pca_wmeta2.to_csv(f"{output_folder}/pca_wmeta_{suffix}.csv", index=False)

pca_corr2 = pca_wmeta2[
    [
        "meer_metage",
        "pet_metage",
        *["PC" + str(i) for i in range(1, pca_result2.shape[1] + 1)],
    ]
].corr()
pca_corr2.to_csv(f"{output_folder}/pca_corr_{suffix}.csv")

metage_explained_var2 = (
    pca_wmeta2[["meer_metage", "pet_metage"]].var() / pca2.explained_variance_.sum()
)
metage_explained_var2.to_csv(
    f"{output_folder}/metage_explained_variance_ratio_{suffix}.csv"
)


weights = meer_clock.set_index("index")["Weight"]
corrs = []
for i, pc in enumerate(pca2.components_):
    pc = pd.Series(index=pca_data2.columns, data=pc)
    shared_weights = pc.index[pc.index.isin(weights.index)]
    corrs.append(
        pd.Series(
            index=["pearsonr", "pval", "spearmanr", "pval", "shared_weights"],
            data=[
                *pearsonr(weights.loc[shared_weights], pc.loc[shared_weights]),
                *spearmanr(weights.loc[shared_weights], pc.loc[shared_weights]),
                len(shared_weights),
            ],
            name=f"PC{i+1}",
        )
    )
pd.concat(corrs, axis=1).T.to_csv(f"{output_folder}/meer_pca_weights_corr_{suffix}.csv")

weights = pet_clock.set_index("index")["Weight"]
corrs = []
for i, pc in enumerate(pca2.components_):
    pc = pd.Series(index=pca_data2.columns, data=pc)
    shared_weights = pc.index[pc.index.isin(weights.index)]
    corrs.append(
        pd.Series(
            index=["pearsonr", "pval", "spearmanr", "pval", "shared_weights"],
            data=[
                *pearsonr(weights.loc[shared_weights], pc.loc[shared_weights]),
                *spearmanr(weights.loc[shared_weights], pc.loc[shared_weights]),
                len(shared_weights),
            ],
            name=f"PC{i+1}",
        )
    )
pd.concat(corrs, axis=1).T.to_csv(f"{output_folder}/pet_pca_weights_corr_{suffix}.csv")


fig, ax = plt.subplots(figsize=(10, 10))
fig = sns.scatterplot(x="PC1", y="meer_metage", data=pca_wmeta2, ax=ax)
plt.savefig(
    f"{output_folder}/pc1_meermetage_{suffix}.png",
    dpi=300,
    bbox_inches="tight",
)

fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x="PC1", y="pet_metage", data=pca_wmeta2, ax=ax)
plt.savefig(
    f"{output_folder}/pc1_petmetage_{suffix}.png",
    dpi=300,
    bbox_inches="tight",
)

fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x="PC1", y="PC2", hue="age_days", data=pca_wmeta2, ax=ax)
plt.savefig(
    f"{output_folder}/pc1_pc2_{suffix}.png",
    dpi=300,
    bbox_inches="tight",
)


fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x="PC1", y="PC2", hue="strain/condition", data=pca_wmeta2, ax=ax)
plt.savefig(
    f"{output_folder}/pc1_pc2_{suffix}_TEST.png",
    dpi=300,
    bbox_inches="tight",
)
