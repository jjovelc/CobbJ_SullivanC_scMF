#!/bin/bash

source ~/activate_my_env.sh


# Define variables
ID="mouse_forelimb"
FASTQ_DIR="/work/vetmed_data/jj/projects/carlySullivan/fastq_files/cellRanger"
REFERENCE="/work/vetmed_data/jj/projects/carlySullivan/fastq_files/cellRanger/mm10"

# Check if output directory already exists
if [ -d "$ID" ]; then
    echo "Error: Output directory '$ID' already exists. Remove or rename before re-running."
    exit 1
fi

# Array of sample names (without _R1/_R2)
SAMPLES=("10_5_cKO" "10_5_WT" "11_5_cKO" "11_5_WT")

# Run CellRanger count for each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"
    
    # Run CellRanger Count
    cellranger count \
        --id="${SAMPLE}_counts" \
        --fastqs="$FASTQ_DIR" \
        --sample="$SAMPLE" \
        --transcriptome="$REFERENCE" \
        --create-bam true

    echo "Completed CellRanger count for sample: $SAMPLE"
done

# Final message
echo "All samples have been processed successfully!"
