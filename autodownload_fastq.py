import subprocess

# # read SRR_Acc_List.txt file with SRA numbers
# with open("SRR_Acc_List.txt", "r") as f:
#     sra_numbers = [line.strip() for line in f]
# # this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
# for sra_id in sra_numbers:
#     print ("Currently downloading: " + sra_id)
#     prefetch = "prefetch " + sra_id
#     print ("The command used was: " + prefetch)
#     subprocess.call(prefetch, shell=True)


with open("SRR_Acc_List-controlsubset.txt", "r") as f:
    sra_numbers = [line.strip() for line in f]
# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print("Generating fastq for: " + sra_id)
    fastq_dump = (
        "fastq-dump --outdir ./data/raw_data --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ./"
        + sra_id
        + "/"
        + sra_id
        + ".sra"
    )
    print("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
