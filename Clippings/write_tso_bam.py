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

def write_degraded_bam_multi(
    BAM, 
    OUTBAM_FILENAME, 
    threads: int = 6):
    '''
    Writes a subsetted BAM file containing only reads with the ts:i tag
    '''
    from multiprocessing import Process
    from multiprocessing import cpu_count

    #print(threads, type(threads))
    threads = int(threads)
    num_cores = int(cpu_count())
    if threads > num_cores:
        threads = num_cores-1
    sys.stdout.write(f'Using {threads} threads. \n')
    
    if os.path.exists(OUTBAM_FILENAME):
        BAD_OUTBAM_FILENAME = OUTBAM_FILENAME
        OUTBAM_FILENAME = OUTBAM_FILENAME.split('.bam')[0] + '-1.bam'
        sys.stderr.write(f'Output BAM {BAD_OUTBAM_FILENAME} exists. using {OUTBAM_FILENAME} \n')
        
    def filesize(file):
        return np.round(os.path.getsize(file) / 1024**2,2)
    
    def chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]
    
    def fetchTSOs(chrom,BAM):
        bam = pysam.AlignmentFile(BAM,'rb')
        TMPBAMFILE = os.path.realpath(str(chrom) + "_tmp.bam")
        TMPBAM = pysam.AlignmentFile(TMPBAMFILE, "wb", template=bam)
        
        Itr = bam.fetch(
            chrom,
            multiple_iterators=True, 
            until_eof=True,
        )
        
        for read in Itr:
            if read.has_tag("ts:i"):
                TMPBAM.write(read)                
        TMPBAM.close()
        # Write a logging function soon
        #sys.stderr.write(f'{os.path.basename(TMPBAMFILE)}: {filesize(TMPBAMFILE)} MB \n')
        
    alignments = pysam.AlignmentFile(BAM, "rb")
    chroms = alignments.header.references
    
    jobs = [Process(target=fetchTSOs(chrom, BAM)) for chrom in chroms]
    for i in chunks(jobs,threads):
        for j in i:
            j.start()

    tmpfiles = [str(chrom) + "_tmp.bam" for chrom in chroms]
    pysam.cat('-o', OUTBAM_FILENAME, *tmpfiles)
    
    for t in tmpfiles:
        os.remove(t)
    
    # Write BAM index
    if os.path.exists(OUTBAM_FILENAME):
        sys.stdout.write(f'Generating BAM index... for {OUTBAM_FILENAME} \n')
        pysam.index(OUTBAM_FILENAME)
        sys.stdout.write(f'Output file: {OUTBAM_FILENAME} \n')
        file_size = filesize(OUTBAM_FILENAME)
        sys.stdout.write(f'File size = {file_size} MB\n')
        
    alignments.close()

def main(args):
    outfile = args.outfile
    BAM = args.BAMFILE
    threads = args.threads
    write_degraded_bam_multi(BAM, outfile, threads=threads)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('BAMFILE', \
        help='Input bam file')
    parser.add_argument('--o', dest='outfile', \
        help='Output filename', default='TSO_reads.bam')
    parser.add_argument('--threads', dest='threads', \
        help='Number of CPU threads to use to speed up processing. Default=4', default=4)
    args = parser.parse_args(sys.argv[1:])
    main(args)
