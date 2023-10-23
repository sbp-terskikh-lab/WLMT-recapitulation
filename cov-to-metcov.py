# This script converts files produced by Bismark v0.15.0 to a format with following columns:
# 1. Chormosome name
# 2. Unique CpG site index which is based on chromosome number and the position on this chromosome
# 3. Coverage
# 4. Methylation level

# Example:
# chr18,chr18_3101849,133,83.4586466165414

import pandas as pd
import os
import subprocess
import concurrent.futures

# results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10122023"
results_folder = "/cjc/miel/martin/WLMT-recapitulation/data/results/results_10102023"

# locations_of_samples_of_interest is a table where the first row contains name of a sample
# and the second - path to the corresponding *.bismark.cov.gz file produced by Bismark v0.15.0

names = pd.read_csv(
    os.path.join(results_folder, "cov_samples_of_interest copy.csv"),
    header=None,
    sep=",",
)

os.makedirs(os.path.join(results_folder, "metcov"), exist_ok=True)


def process_sample(name, file):
    print("Starting processing ", name)
    if not file.endswith(".gz"):
        raise ValueError("File is not gzipped")
    else:
        data0 = pd.read_csv(
            os.path.join(results_folder, "methylation", file),
            compression="gzip",
            header=None,
            sep="\t",
        )
    index = data0[0].astype(str) + "_" + data0[1].astype(str)
    coverage = data0[4] + data0[5]
    data1 = pd.DataFrame(
        {
            "chr": data0[0].astype(str),
            "index": index,
            "coverage": coverage,
            "metlev": data0[3],
        }
    )
    outfilename = os.path.join(results_folder, "metcov", f"{name}.metcov.gz")
    data1.to_csv(outfilename, compression="gzip", index=False)
    return "Finished processing " + name + " and saved to " + outfilename + "\n


with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = [
        executor.submit(process_sample, name, file) for name, file in names.values
    ]
    for future in futures:
        print(future.result())
    executor.shutdown(wait=True)
