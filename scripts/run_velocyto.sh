#!/usr/bin/bash


DIR=$1 
# e.g. 10_5_cKO_counts (cellranger count folder) 

GENES=$2
# or GENES="/work/vetmed_data/jj/projects/carlySullivan/fastq_files/cellRanger/mm10/genes/genes.gtf"

velocyto run10x -@ 8 "$DIR" "$GENES"
