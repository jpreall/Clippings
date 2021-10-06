"""
write_tso_bam.py

Writes a subsetted BAM file containing only reads with the ts:i tag

Inputs:
    BAMFILE = BAM file produced by Cellranger v4+ with TSO reads tagged with ts:i
	optional: 
    	--out = specific desired output folder.  defaults to same is input file
    Usage: Clippings_count_miRNA.py BAMFILE

Author: jpreall@cshl.edu
Date: 2021-08
"""

import os, sys
import argparse
import pysam
import numpy as np


def write_degraded_bam(full_path_to_BAM, OUTBAM_FILENAME):
    
	'''Writes a subsetted BAM file containing only reads with the ts:i tag
	There is probably a good way to parallelize this in the future
	'''
    
	BAM = os.path.basename(full_path_to_BAM)  
	print('Extracting TSO-containing reads from',BAM,'...')
    
	alignments = pysam.AlignmentFile(BAM, "rb")
	OUTBAM = pysam.AlignmentFile(OUTBAM_FILENAME, "wb", template=alignments)
    
	tally = 0
    
	for read in alignments.fetch(until_eof=True):
		tally += 1
		if tally % 1e7 == 0:
			print('Lines read:',f'{tally:,}')
			
		if read.has_tag("ts:i"):
			OUTBAM.write(read)

	OUTBAM.close()
	alignments.close()


    ## Write the .bai index file
	if os.path.exists(OUTBAM_FILENAME):
		print('Generating BAM index...')
		pysam.index(OUTBAM_FILENAME)
		
		print('Output file:',OUTBAM_FILENAME)
		print('File size = ',np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2,2),'MB')	
		
	OUTBAM.close()
        
        
def main(args):
	outpath = args.out ## NEED A FRIENDLY WAY TO CHECK IF PATH EXISTS
	BAM = args.BAMFILE
	full_path_to_BAM = os.path.realpath(BAM)
	
	if outpath == None:
		outpath = os.path.dirname(full_path_to_BAM)
		
	if not os.path.isdir(outpath):
		print('Warning: output path not found')
        
	OUTBAM_FILENAME = os.path.join(outpath,'TSO_reads.bam')
	
	write_degraded_bam(full_path_to_BAM, OUTBAM_FILENAME)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('--out', dest='out', help='Output folder', default=None)
    

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    main(args)