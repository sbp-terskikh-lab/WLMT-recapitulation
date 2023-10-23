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


# Initialize WD
root_folder = "/cjc/miel/martin/WLMT-recapitulation"
results_folder = f"{root_folder}/data/results/results_10102023"
# results_folder = f"{root_folder}/data/results/results_10122023"
output_folder = f"{results_folder}/analysis_result"
os.makedirs(output_folder, exist_ok=True)

# # Read in metadata
# petkovich_meta = pd.read_excel(
#     f"{root_folder}/sample_meta_NIHMS863021-supplement-2.xlsx", header=3
# )
# gsm_meta = pd.read_csv(f"{root_folder}/metamatch.txt", header=None)
# gsm_meta.columns = ["Sample Name", "Sample"]
# sra_table = pd.read_csv(f"{root_folder}/SraRunTable-2.txt")
# nih_meta = pd.merge(gsm_meta, sra_table, on="Sample Name", how="inner").rename(
#     columns={"Sample": "Sample Name", "Sample Name": "GEO_Accession"}
# )
# all_meta = pd.merge(nih_meta, petkovich_meta, on="Sample Name", how="inner")
# all_meta.columns = all_meta.columns.str.lower()
# all_meta["age_days"] = all_meta["age"].astype(float) * 30
# all_meta.to_csv(f"{output_folder}/all_meta.csv", index=False)

all_meta = pd.read_csv(f"{output_folder}/all_meta.csv")


# # Read in Data
# files = os.listdir(f"{results_folder}/metlev")
# files = [f for f in files if f.endswith("full_autosome.csv.gz")]

# files = files[:20]

# def read_csv(file):
#     return pd.read_csv(f"{results_folder}/metlev/{file}")

# metlev_data = []
# with ThreadPoolExecutor() as executor:
#     metlev_data = list(executor.map(read_csv, files))


# Read in WLMT clock
xlsx = pd.ExcelFile(f"{root_folder}/clock_weights_elife-40675-supp3-v2.xlsx")
meer_clock = pd.read_excel(
    xlsx, sheet_name="Whole lifespan multi-tissue", header=0
).dropna(axis=0, how="any")
meer_clock["chr"] = meer_clock["Chromosome"].str.replace("chr", "").astype(int)
meer_clock.rename(columns={"Position": "position"}, inplace=True)
meer_clock["index"] = (
    meer_clock[["chr", "position"]]
    .astype(int)
    .astype(str)
    .apply(lambda x: "_".join(x), axis=1)
)


def calculate_metage(i):
    mdat = metlev_data[i]
    shared = np.intersect1d(meer_clock["index"].values, mdat["index"].values)
    relevant_data = mdat.set_index("index").loc[shared]
    metage = relevant_data.T.dot(
        meer_clock.set_index("index").loc[shared, "Weight"]
    ) / np.linalg.norm(meer_clock["Weight"])
    metage = metage.to_frame().reset_index().rename(columns={0: "metage"})
    name = metage["index"].str.split("_", expand=True)[0]
    metage["run"] = name
    return name, metage, shared, len(shared)


# names = []
# metage = []
# shared = []
# num_overlap = []
# with ProcessPoolExecutor() as executor:
#     futures = [
#         executor.submit(calculate_metage, i) for i in np.arange(len(metlev_data))
#     ]
#     for future in as_completed(futures):
#         name, metage_, shared_, num_overlap_ = future.result()
#         names.append(name)
#         metage.append(metage_)
#         shared.append(shared_)
#         num_overlap.append(num_overlap_)
# metage = pd.concat(metage, ignore_index=True)
# metage.to_csv(f"{output_folder}/meer_metage.csv", index=False)
# np.save(f"{output_folder}/meer_shared.npy", np.array(shared))
# np.save(f"{output_folder}/meer_num_overlap.npy", np.array(num_overlap))
# np.save(f"{output_folder}/meer_names.npy", np.array(names))

metage = pd.read_csv(f"{output_folder}/meer_metage.csv")

metage.rename(columns={"metage": "meer_metage"}, inplace=True)
metage_wmeta = pd.merge(all_meta, metage, on="run", how="inner")


# fig, ax = plt.subplots(figsize=(10, 10))
# sns.scatterplot(x="age_days", y="meer_metage", data=metage_wmeta, ax=ax)
# plt.savefig(
#     f"{output_folder}/age_meermetage.png",
#     dpi=300,
#     bbox_inches="tight",
# )
# agecorr = pearsonr(metage_wmeta["age_days"], metage_wmeta["meer_metage"])
# np.save(
#     f"{output_folder}/meer_regstats.npy", {"pearsonr": agecorr[0], "pval": agecorr[1]}
# )


pet_clock = pd.read_excel(
    f"{root_folder}/clock_weights_NIHMS863021-supplement-4.xlsx", header=3
)
pet_clock["chr"] = pet_clock["Chromosome"].str.replace("chr", "").astype(int)
pet_clock.rename(columns={"Position": "position"}, inplace=True)
pet_clock["index"] = (
    pet_clock[["chr", "position"]]
    .astype(int)
    .astype(str)
    .apply(lambda x: "_".join(x), axis=1)
)
pet_clock


# def calculate_metage(i):
#     mdat = metlev_data[i]
#     shared = np.intersect1d(pet_clock["index"].values, mdat["index"].values)
#     relevant_data = mdat.set_index("index").loc[shared]
#     metage = relevant_data.T.dot(
#         pet_clock.set_index("index").loc[shared, "Weight"]
#     ) / np.linalg.norm(pet_clock["Weight"])
#     metage = metage.to_frame().reset_index().rename(columns={0: "metage"})
#     name = metage["index"].str.split("_", expand=True)[0]
#     metage["run"] = name
#     return name, metage, shared, len(shared)


# names = []
# metage = []
# shared = []
# num_overlap = []
# with ProcessPoolExecutor() as executor:
#     futures = [
#         executor.submit(calculate_metage, i) for i in np.arange(len(metlev_data))
#     ]
#     for future in as_completed(futures):
#         name, metage_, shared_, num_overlap_ = future.result()
#         names.append(name)
#         metage.append(metage_)
#         shared.append(shared_)
#         num_overlap.append(num_overlap_)
# metage = pd.concat(metage, ignore_index=True)
# metage.to_csv(f"{output_folder}/pet_metage.csv", index=False)
# np.save(f"{output_folder}/pet_shared.npy", np.array(shared))
# np.save(f"{output_folder}/pet_num_overlap.npy", np.array(num_overlap))
# np.save(f"{output_folder}/pet_names.npy", np.array(names))

metage = pd.read_csv(f"{output_folder}/pet_metage.csv")

metage.rename(columns={"metage": "pet_metage"}, inplace=True)
metage_wmeta = pd.merge(metage_wmeta, metage, on="run", how="inner")

# fig, ax = plt.subplots(figsize=(10, 10))
# sns.scatterplot(x="age_days", y="pet_metage", data=metage_wmeta, ax=ax)
# plt.savefig(
#     f"{output_folder}/age_petmetage.png",
#     dpi=300,
#     bbox_inches="tight",
# )
# agecorr = pearsonr(metage_wmeta["age_days"], metage_wmeta["pet_metage"])
# np.save(
#     f"{output_folder}/pet_regstats.npy", {"pearsonr": agecorr[0], "pval": agecorr[1]}
# )


# full_idx = list(set().union(*[set(df["index"].values) for df in metlev_data]))
# genelistmaxidx = len(full_idx)
# gene_to_idx = dict(zip(full_idx, np.arange(len(full_idx))))
# with open(f"{output_folder}/gene_to_idx.pickle", "wb") as f:
#     pickle.dump(gene_to_idx, f, protocol=pickle.HIGHEST_PROTOCOL)

with open(f"{output_folder}/gene_to_idx.pickle", "rb") as f:
    gene_to_idx = pickle.load(f)
    full_idx = list(gene_to_idx.keys())
    genelistmaxidx = len(full_idx)


# def metlev_to_sparse(df):
#     # Get name, unpack data
#     data_colname = df.columns[1]
#     name = data_colname.strip("_metlev")
#     genes = df["index"].values
#     geneidxs = df["index"].map(gene_to_idx).values
#     data = df[data_colname].values
#     df = None
#     # Create sparse matrix
#     spm = coo_matrix(
#         (data, (geneidxs, np.zeros(len(geneidxs), dtype=int))),
#         shape=(genelistmaxidx, 1),
#     )
#     data = None
#     geneidxs = None
#     return name, genes, spm


# names, sample_genes, matrices = [], {}, []
# # with ThreadPoolExecutor(max_workers=20) as executor:
# with ProcessPoolExecutor(max_workers=100) as executor:
#     futures = [executor.submit(metlev_to_sparse, df) for df in metlev_data]
#     for future in as_completed(futures):
#         name, genes, spm = future.result()
#         names.append(name)
#         sample_genes[name] = genes
#         matrices.append(spm)
#         futures.remove(future)
#         print(f"Finished {name}")
#         print(f"Remaining: {len(futures)}")
#     executor.shutdown(wait=True)

# metlev_merge = hstack(matrices)
# metlev_merge = pd.DataFrame(
#     data=metlev_merge.toarray(),
#     index=full_idx,
#     columns=names,
# )
# metlev_merge.to_hdf(f"{output_folder}/metlev_merge.hdf5", key="df", format="table")
# joblib.dump(sample_genes, f"{output_folder}/sample_genes.joblib")


metlev_merge = pd.read_hdf(f"{output_folder}/metlev_merge.hdf5", key="df")
sample_genes = joblib.load(f"{output_folder}/sample_genes.joblib")

# # Open an HDF5 file in write mode
# with tables.open_file(f'{output_folder}/metlev_combined.h5', mode='w') as f:
#     # Create a group in the file
#     group = f.create_group('/', 'mygroup')

#     # Create a dataset in the group
#     dset = f.create_carray(group, 'mydata', obj=metlev_merge.toarray())
# np.save(f'{output_folder}/metlev_combined_names.npy', np.array(names))


# # Open the HDF5 file in read mode
# with tables.open_file(f"{output_folder}/metlev_combined.h5", mode="r") as f:
#     # Get the dataset from the file
#     dset = f.get_node("/mygroup/mydata")

#     # Read the data from the dataset
#     metlev_merge = dset.read()
# names = np.load(f"{output_folder}/metlev_combined_names.npy", allow_pickle=True)

# metlev_merge = pd.DataFrame(
#     metlev_merge, columns=names,
# )

mean = metlev_merge.mean(axis=1)
pca_data = metlev_merge.sub(mean, axis=0)
pca_data = pca_data.T
# pca = PCA(n_components=pca_data.shape[0])
# pca_result = pca.fit_transform(pca_data)
# joblib.dump(pca, f"{output_folder}/pca.joblib")
# # with open(f"{output_folder}/pca.pickle", "wb") as f:
# #     pickle.dump(pca, f)
# np.save(
#     f"{output_folder}/pca_explained_variance_ratio.npy", pca.explained_variance_ratio_
# )

# pca_result = pd.DataFrame(
#     pca_result,
#     index=pca_data.index,
#     columns=["PC" + str(i) for i in range(1, pca_result.shape[0] + 1)],
# )
# pca_wmeta = pd.merge(
#     metage_wmeta, pca_result, left_on="run", right_index=True, how="inner"
# )
# pca_wmeta.to_csv(f"{output_folder}/pca_wmeta.csv", index=False)

# pca_corr = pca_wmeta[
#     [
#         "meer_metage",
#         "pet_metage",
#         *["PC" + str(i) for i in range(1, pca_result.shape[1] + 1)],
#     ]
# ].corr()
# pca_corr.to_csv(f"{output_folder}/pca_corr.csv")

# metage_explained_var = (
#     pca_wmeta[["meer_metage", "pet_metage"]].var() / pca.explained_variance_.sum()
# )
# metage_explained_var.to_csv(f"{output_folder}/metage_explained_variance_ratio.csv")


# weights = meer_clock.set_index("index")["Weight"]
# corrs = []
# for i, pc in enumerate(pca.components_):
#     pc = pd.Series(index=pca_data.columns, data=pc)
#     shared_weights = pc.index[pc.index.isin(weights.index)]
#     corrs.append(
#         pd.Series(
#             index=["pearsonr", "pval", "spearmanr", "pval", "shared_weights"],
#             data=[
#                 *pearsonr(weights.loc[shared_weights], pc.loc[shared_weights]),
#                 *spearmanr(weights.loc[shared_weights], pc.loc[shared_weights]),
#                 len(shared_weights),
#             ],
#             name=f"PC{i+1}",
#         )
#     )
# pd.concat(corrs, axis=1).T.to_csv(f"{output_folder}/meer_pca_weights_corr.csv")

# weights = pet_clock.set_index("index")["Weight"]
# corrs = []
# for i, pc in enumerate(pca.components_):
#     pc = pd.Series(index=pca_data.columns, data=pc)
#     shared_weights = pc.index[pc.index.isin(weights.index)]
#     corrs.append(
#         pd.Series(
#             index=["pearsonr", "pval", "spearmanr", "pval", "shared_weights"],
#             data=[
#                 *pearsonr(weights.loc[shared_weights], pc.loc[shared_weights]),
#                 *spearmanr(weights.loc[shared_weights], pc.loc[shared_weights]),
#                 len(shared_weights),
#             ],
#             name=f"PC{i+1}",
#         )
#     )
# pd.concat(corrs, axis=1).T.to_csv(f"{output_folder}/pet_pca_weights_corr.csv")


# fig, ax = plt.subplots(figsize=(10, 10))
# fig = sns.scatterplot(x="PC1", y="meer_metage", data=pca_wmeta, ax=ax)
# plt.savefig(
#     f"{output_folder}/pc1_meermetage.png",
#     dpi=300,
#     bbox_inches="tight",
# )

# fig, ax = plt.subplots(figsize=(10, 10))
# sns.scatterplot(x="PC1", y="pet_metage", data=pca_wmeta, ax=ax)
# plt.savefig(
#     f"{output_folder}/pc1_petmetage.png",
#     dpi=300,
#     bbox_inches="tight",
# )

# fig, ax = plt.subplots(figsize=(10, 10))
# sns.scatterplot(x="PC1", y="PC2", hue="age_days", data=pca_wmeta, ax=ax)
# plt.savefig(
#     f"{output_folder}/pc1_pc2.png",
#     dpi=300,
#     bbox_inches="tight",
# )

propnan_thresh = 0.01
suffix = f"{propnan_thresh}propnan_filtered_16sample"
prop_nan = (metlev_merge == 0).sum(axis=1) / metlev_merge.shape[1]
pca_data2_ = pca_data.loc[:, prop_nan <= propnan_thresh]
pca_data2_.loc[:, meer_clock["index"]] = pca_data.loc[:, meer_clock["index"]]
np.savetxt(
    f"{output_folder}/num_features_{suffix}.csv", pca_data2_.shape, delimiter=","
)
prop_nan.rename("prop_nan")
prop_nan.reset_index().to_csv(f"{output_folder}/prop_nan_{suffix}.csv", index=False)
output_folder = f"{output_folder}/pca_{suffix}"
os.makedirs(output_folder, exist_ok=True)


def run_pca2(pca_data2, boot):
    suffix = f"{propnan_thresh}propnan_filtered_16sample_boot{boot}"
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

    # with open(f"{output_folder}/pca_{suffix}.pickle", "wb") as f:
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
    pd.concat(corrs, axis=1).T.to_csv(
        f"{output_folder}/meer_pca_weights_corr_{suffix}.csv"
    )

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
    pd.concat(corrs, axis=1).T.to_csv(
        f"{output_folder}/pet_pca_weights_corr_{suffix}.csv"
    )

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


with ProcessPoolExecutor(max_workers=20) as executor:
    futures = [
        executor.submit(run_pca2, pca_data2_.sample(16, random_state=boot), boot)
        for boot in range(100)
    ]
    for future in tqdm(as_completed(futures)):
        future.result()
    executor.shutdown(wait=True)
