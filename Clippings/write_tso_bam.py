"""
write_tso_bam.py

Writes a subsetted BAM file containing only reads with the ts:i tag

Inputs:
    BAMFILE = BAM file produced by Cellranger v4+ with TSO reads tagged with ts:i
    optional:
        --o = output filename.  Default: 'TSO_reads.bam'
    Usage: python write_tso_bam.py BAMFILE
Author: bhe@cshl.edu
Date: 2021-08
"""

import os
import sys
import argparse
import pysam
import numpy as np

def get_template(full_path_to_BAM):
    BAM = full_path_to_BAM
    alignments_template = pysam.AlignmentFile(BAM, "rb")
    return alignments_template

def write_degraded_bam_simple(BAM, OUTBAM_FILENAME):
    '''
    Writes a subsetted BAM file containing only reads with the ts:i tag
	There is probably a good way to parallelize this in the future
	'''
    
    sys.stdout.write(f'Extracting TSO-containing reads from {BAM} ... \n')
    alignments = pysam.AlignmentFile(BAM, "rb")
    OUTBAM = pysam.AlignmentFile(OUTBAM_FILENAME, "wb", template=alignments)
    
    tally = 0
    
    for read in alignments.fetch(until_eof=True):
        tally += 1
        if tally % 1e7 == 0:
            sys.stdout.write(f'Lines read: {tally:,} \n')
			
        if read.has_tag("ts:i"):
            OUTBAM.write(read)

    OUTBAM.close()
    alignments.close()

    if os.path.exists(OUTBAM_FILENAME):
        sys.stdout.write('Generating BAM index... \n')
        pysam.index(OUTBAM_FILENAME)
        sys.stdout.write(f'Output file: {OUTBAM_FILENAME} \n')
        file_size = np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2,2)
        sys.stdout.write(f'File size = {file_size} MB \n')

    OUTBAM.close()

def main(args):
    outfile = args.outfile
    BAM = args.BAMFILE
    write_degraded_bam_simple(BAM, outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('--o', dest='outfile', help='Output filename', default='TSO_reads.bam')
    args = parser.parse_args(sys.argv[1:])
    main(args)
