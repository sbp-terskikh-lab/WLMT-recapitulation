#!/bin/bash

# Set the input variables
N_CORES=12
N_PARALLEL=8
DATA_FOLDER="/cjc/miel/martin/WLMT-recapitulation/data"
RESULTS_FOLDER="$DATA_FOLDER/results/results_10102023"
TEMP_FOLDER="$RESULTS_FOLDER/temp"
LOCK_FOLDER="$RESULTS_FOLDER/lockfiles"
SAMPLE_SHEET="$DATA_FOLDER/raw_data/SampleSheet.csv"
GENOME_FOLDER="$DATA_FOLDER/genome_folder/grcm38"
RAW_DATA_DIR="$DATA_FOLDER/raw_data"
TRIM_OUTPUT_DIR="$RESULTS_FOLDER/trimmed"
ALIGN_OUTPUT_DIR="$RESULTS_FOLDER/alignments"
OUTPUT_DIR3="$RESULTS_FOLDER/methylation"
LOCK_FILE="$LOCK_FOLDER/lock"

# Generate a unique ID for this script instance
SCRIPT_ID=$$

# Randomize the sample sheet list
tail -n +2 "$SAMPLE_SHEET" | shuf > "$TEMP_FOLDER/sample_sheet_shuf.$SCRIPT_ID.csv"

# Process each sample
while IFS=, read -r sample fastq_1 fastq_2; do
  echo "Processing sample $sample: $fastq_1, $fastq_2"

  # Check if the lock file exists for this sample
  if [ -f "$LOCK_FILE.$sample" ]; then
    echo "Skipping sample $sample (already running or completed)"
    continue
  fi

  # Create the lock file for this sample
  touch "$LOCK_FILE.$sample"
  
  # Run trim_galore
  if ! trim_galore --rrbs --paired --gzip --three_prime_clip_R1 15 --three_prime_clip_R2 15 --cores $N_CORES -o "$TRIM_OUTPUT_DIR" "$RAW_DATA_DIR/$fastq_1" "$RAW_DATA_DIR/$fastq_2"; then
    echo "Error processing sample $sample: trim_galore failed"
    rm -f "$LOCK_FILE.$sample"
    continue
  fi

  # Run bismark
  if ! bismark --parallel $N_PARALLEL --multicore 4 --temp_dir "$ALIGN_OUTPUT_DIR/temp" -o "$ALIGN_OUTPUT_DIR" --gzip -genome "$GENOME_FOLDER" -1 "$TRIM_OUTPUT_DIR/${fastq_1%.fastq.gz}_val_1.fq.gz" -2 "$TRIM_OUTPUT_DIR/${fastq_2%.fastq.gz}_val_2.fq.gz"; then
    echo "Error processing sample $sample: bismark failed"
    rm -f "$TRIM_OUTPUT_DIR/${fastq_1%.fastq.gz}_val_1.fq.gz" "$TRIM_OUTPUT_DIR/${fastq_2%.fastq.gz}_val_2.fq.gz"
    rm -f "$LOCK_FILE.$sample"
    continue
  fi

  # Remove intermediate files
  rm -f "$TRIM_OUTPUT_DIR/${fastq_1%.fastq.gz}_val_1.fq.gz" "$TRIM_OUTPUT_DIR/${fastq_2%.fastq.gz}_val_2.fq.gz"

  # Run bismark_methylation_extractor
  if ! bismark_methylation_extractor -p --comprehensive --bedGraph --merge_non_cpg --gzip --multicore $N_CORES -o "$OUTPUT_DIR3" "$ALIGN_OUTPUT_DIR/${fastq_1%.fastq.gz}_val_1_bismark_bt2_pe.bam"; then
    echo "Error processing sample $sample: bismark_methylation_extractor failed"
    rm -f "$ALIGN_OUTPUT_DIR/${fastq_1%.fastq.gz}_val_1_bismark_bt2_pe.bam"
    rm -f "$LOCK_FILE.$sample"
    continue
  fi

  # Remove intermediate files
  rm -f "$ALIGN_OUTPUT_DIR/${fastq_1%.fastq.gz}_val_1_bismark_bt2_pe.bam"
  
done < "$TEMP_FOLDER/sample_sheet_shuf.$SCRIPT_ID.csv"

# Remove the shuffled sample sheet file
rm -f "$TEMP_FOLDER/sample_sheet_shuf.$SCRIPT_ID.csv"