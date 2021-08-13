# Clippings
Use pysam to walk through the bam file and count every read with a ts:i tag, indicating the tso was soft clipped.
Count mapped reads harboring TSO sub-sequene.

## Overview

## Installation
### 1. Create a new directory and clone repository
```
(your_env) $ git clone https://github.com/jpreall/Clippings.git
(your_env) $ pip install -e Clippings

## Example Usage
(your_env) $ qsub -R y -l m_mem_free=8G -pe threads 16 deg_count_with_UMIs.py [path_to_bam_file] --outfile [output_dictionary_textName]  --cores 12 --TSSgtf [path_to_gtf] --out [output_folder_name] --genome [path_to_genome]


