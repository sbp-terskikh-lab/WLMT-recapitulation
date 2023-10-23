import os
import pandas as pd

raw_data_files = pd.Series(os.listdir("./data/raw_data"))
raw_data_files = raw_data_files[raw_data_files.str.contains(".fastq.gz")]
sra_numbers = raw_data_files.str.extract(r"(SRR\d+)", expand=False)
pass_number = raw_data_files.str.extract(r"(pass_\d+)", expand=False)
samplesheet = (
    pd.DataFrame(
        {
            "sra_number": sra_numbers,
            "pass_number": pass_number,
            "filepath": raw_data_files,
        }
    )
    .pivot(index="sra_number", columns="pass_number", values="filepath")
    .reset_index()
)
samplesheet.columns = ["sample", "fastq_1", "fastq_2"]
samplesheet.to_csv("SampleSheet.csv", index=False)
