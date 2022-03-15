"""
write_tso_bam.py

Writes a subsetted BAM file containing only reads with the ts:i tag

Inputs:
    BAMFILE = BAM file produced by Cellranger v4+ with TSO reads tagged with ts:i
    optional:
        --out = specific desired output folder.  defaults to same is input file
    Usage: python write_tso_bam.py BAMFILE
    TODO: doesn't work unless in directory of bam (?) Need to update code
Author: jpreall@cshl.edu
Date: 2021-08
"""

import os
import sys
import argparse
import pysam
import numpy as np
from pararead import ParaReadProcessor
import time
import glob


class TSOReadCounter(ParaReadProcessor):
    """ Read counter for TSO labelled reads. """

    def __call__(self, chromosome, _=None):
        import os
        # doesn't work because count can't handle multiple iterators
        # n_reads = self.readsfile.count(chromosome)

        # Here we want to go through each read and do what we need.
        tso_reads = 0
        reads = self.fetch_chunk(chromosome)

        print('Time started for TSO BAM ({0}) write-out:'.format(chromosome), time.asctime())
        # tso_bamFolder = 'tso_bamFolder'
        # tso_bamFolderPath = os.path.join(os.getcwd(), tso_bamFolder)
        tso_bamFolderPath = os.getcwd()
        TSOBAM_FILENAME = os.path.join(tso_bamFolderPath, str(chromosome) + '_tmp_tso.bam')

        TSOBAM = pysam.AlignmentFile(TSOBAM_FILENAME, "wb", template=get_template(args.BAMFILE))

        for read in reads:
            if read.has_tag("ts:i"):
                TSOBAM.write(read)
                tso_reads += 1
        TSOBAM.close()

        # _LOGGER.debug("Chromosome: '{}'; n_reads {}".format(chromosome, n_reads))
        with open(self._tempf(chromosome), 'w') as f:
            f.write("{}\t{}".format(chromosome, tso_reads))
        return chromosome


def get_template(full_path_to_BAM):
    # BAM = os.path.basename(full_path_to_BAM)
    BAM = full_path_to_BAM
    alignments_template = pysam.AlignmentFile(BAM, "rb")
    return alignments_template


def main(args):
    import time
    print('****MAIN STEP****: START', time.asctime())

    tso_bamFolder = args.out
    import os
    tso_bamFolderPath = os.path.join(os.getcwd(), tso_bamFolder)
    os.mkdir(tso_bamFolderPath)
    os.chdir(tso_bamFolderPath)
    print('tso_bamFolderPath', tso_bamFolderPath)

    # custom outpath currently doesnt work
    # outpath = args.out  # NEED A FRIENDLY WAY TO CHECK IF PATH EXISTS
    BAM = args.BAMFILE
    BAM_align = pysam.AlignmentFile(BAM, 'rb')

    print("Creating counter")
    counter = TSOReadCounter(args.BAMFILE, cores=args.cores, action="CountReads")
    print('****MAIN STEP****: FINISHED THE READCOUNTER', time.asctime())

    print("Registering files")
    counter.register_files()
    print('****MAIN STEP****: FINISHED REGISTERING', time.asctime())

    print("Counting reads: {}".format(args.BAMFILE))
    good_chromosomes = counter.run()
    print('****MAIN STEP****: FINISHED counter.run()', time.asctime())

    # print("Collecting read counts: {}".format(args.outdictionary))
    counter.combine(good_chromosomes, chrom_sep="\n")
    print('****MAIN STEP****: FINISHED COMBINING', time.asctime())

    print('Combine bam files', time.asctime())

    print("Concatenating bam files", time.asctime())
    # Now that the chromosome bams are written out, we can proceed
    allBam = '*.bam'
    bamFiles = glob.glob(os.path.join(os.getcwd(), allBam))
    # print('bamFiles: ', bamFiles)

    # Get the ordered list of chromosomes based on original bam
    chr_list_order = list()
    for chromosome in BAM_align.header['SQ']:
        # print(chromosome['SN'] + '_tmp_tso.bam')
        chr_list_order.append(chromosome['SN'] + '_tmp_tso.bam')

    # next get the full path for each chromosome
    chr_path_list_order = list()
    for i in chr_list_order:
       # print(os.path.join(os.getcwd(), i))
        chr_path_list_order.append(os.path.join(os.getcwd(), i))

    # then get the subset of chromosome bams that exist
    # in the folder from chr_path_list_order
    # https://stackoverflow.com/questions/23529001/ordered-intersection-of-two-lists-in-python
    setOfThingsToKeep = frozenset(bamFiles)
    intersectionKeepOrder = [x for x in chr_path_list_order if x in setOfThingsToKeep]

    # finally concatenate the ts:i labelled bam files together
    pysam.cat('-o', 'TSO_reads_bam.bam', *intersectionKeepOrder)
    print('****MAIN STEP****: FINISHED CONCATENATING', time.asctime())

    # Write the .bai index file
    full_tso_bam = os.path.join(os.getcwd(), 'TSO_reads_bam.bam')
    if os.path.exists(full_tso_bam):
        print('Generating BAM index... ', time.asctime())
        pysam.index(full_tso_bam)
        print('Output file:', full_tso_bam)
        print('File size = ', np.round(os.path.getsize(full_tso_bam) / 1024 ** 2, 2), 'MB')
    else:
        print('Path doesn\'t exist')
    print('Finished script', time.asctime())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('--outfolder', dest='out', help='Output folder', default='tso_bamFolder')
    parser.add_argument('--cores', default=10, help="Number of cores.")
    # testing purposes start
#    sys.argv = ['write_tso_bam.py',
#                '/mnt/grid/scc/data/Preall/Preall_CR01/count/Preall_CR01_H_neg/outs/possorted_genome_bam.bam',
#                '--outfolder', 'tso_bamFolder_CR01_H_neg',
#                '--cores', '10']
#    print(sys.argv[1:])
    # testing purposes end
    args = parser.parse_args(sys.argv[1:])
    main(args)
