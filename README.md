# Clippings
Use pysam to walk through the bam file and count every read with a ts:i tag, indicating the tso was soft clipped.
Count mapped reads harboring TSO sub-sequene.

## Overview

## Installation
### 1. Create a new directory and clone repository
```
(your_env) $ git clone https://github.com/jpreall/Clippings.git
(your_env) $ pip install -e Clippings
```

## Example Usage
### 1. deg_count_with_UMIs.py 
##### This will output barcodes.tsv.gz, dictionary.txt, features.tsv.gz, jsonFolder, matrix.mtx.gz, raw_clipped_features_matrix.h5
```
# Highly recommend user to submit as a job if using HPC in order to prevent over overuse of resources. This should take 5-10 minutes.
(your_env) $ python deg_count_with_UMIs.py [path_to_bam_file] --TSSgtf [path_to_gtf/genes.gtf] --cores 10 --outdir [output_folder_name] --mtx [True/False]
```

### 2. Clippings_count_miRNAs.py
This will output: 
  \[SAMPLE\]_sorted_miRNA_reads.bam
  miRNA_read_details.csv.gz
  miRNA_count_matrix.csv.gz
  TSO_distances_to_Drosha.csv.gz
  feature_bc_matrix_with_miRNAs/
```
(your_env) $ python Clippings_count_miRNAs.py [path_to_bam_file] [path_to_gff3] --outdir [output_folder_name] --genome [path_to_genome] --raw [path_to_raw_matrix]
```

### 3. write_tso_bam.py
This will output TSO_reads.bam and TSO_reads.bam.bai
Supports multiple CPU threads to speed up processing. 
```
(your_env) $ python write_tso_bam.py [path_to_bam_file] --threads [threads]
```
