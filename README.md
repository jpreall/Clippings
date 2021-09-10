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
```
# Highly recommend user to submit as a job if using HPC in order to prevent over overuse of resources. This should take 5-10 minutes.
(your_env) $ python deg_count_with_UMIs.py [path_to_bam_file] --TSS [path_to_gtf] --outdir [output_folder_name] --genome [path_to_genome]
```
### 2. Clippings_count_miRNAs.py
```
(your_env) $ python Clippings_count_miRNAs.py [path_to_bam_file] [path_to_gff3] --outdir [output_folder_name] --genome [path_to_genome] --raw [path_to_raw_matrix]
```
