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
        # doesn't work because count can't handle multiple iterators
        # n_reads = self.readsfile.count(chromosome)

        # Here we want to go through each read and do what we need.
        tso_reads = 0
        reads = self.fetch_chunk(chromosome)

        print('Time started for TSO BAM write-out:', time.asctime())
        # 2021.08.09 Writing out dictionary as json to then re-load and combine
        tso_bamFolder = 'tso_bamFolder'
        import os
        tso_bamFolderPath = os.path.join(os.getcwd(), tso_bamFolder)
        # change back if necessary

        OUTBAM_FILENAME = os.path.join(tso_bamFolderPath, str(chromosome) + '_TSO_reads.bam')
        TSOBAM_FILENAME = os.path.join(tso_bamFolderPath, str(chromosome) + '_tmp_tso.bam')
        #alignments = pysam.AlignmentFile(BAM, "rb")

        TSOBAM = pysam.AlignmentFile(TSOBAM_FILENAME, "wb", template=get_template(args.BAMFILE))

        for read in reads:
            if read.has_tag("ts:i"):
                TSOBAM.write(read)
                tso_reads += 1
        TSOBAM.close()

        #_LOGGER.debug("Chromosome: '{}'; n_reads {}".format(chromosome, n_reads))
        with open(self._tempf(chromosome), 'w') as f:
            f.write("{}\t{}".format(chromosome, tso_reads))
        return chromosome


def get_template(full_path_to_BAM):
    BAM = os.path.basename(full_path_to_BAM)
    alignments_template = pysam.AlignmentFile(BAM, "rb")
    return alignments_template


def write_degraded_bam(full_path_to_BAM, TMPBAM_FILENAME, OUTBAM_FILENAME):
    '''Writes a subsetted BAM file containing only reads with the ts:i tag
    There is probably a good way to parallelize this in the future
    '''

    BAM = os.path.basename(full_path_to_BAM)
    print('Extracting TSO-containing reads from', BAM, '...')

    alignments = pysam.AlignmentFile(BAM, "rb")
    TMPBAM = pysam.AlignmentFile(TMPBAM_FILENAME, "wb", template=alignments)

    tally = 0

    for read in alignments.fetch(until_eof=True):
        tally += 1
        if tally % 1e7 == 0:
            print('Lines read:', f'{tally:,}')

        if read.has_tag("ts:i"):
            TMPBAM.write(read)

    TMPBAM.close()
    alignments.close()

    # Write the .bai index file
    if os.path.exists(TMPBAM_FILENAME):

        print('Sorting output...')
        pysam.sort("-o", OUTBAM_FILENAME, TMPBAM_FILENAME)
        os.remove(TMPBAM_FILENAME)
        if not os.path.exists(TMPBAM_FILENAME):
            print('temp files deleted successfully...')
        print('Generating BAM index...')
        pysam.index(OUTBAM_FILENAME)

        print('Output file:', OUTBAM_FILENAME)
        print('File size = ', np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2, 2), 'MB')

    TMPBAM.close()


def main(args):
    import time
    print('****MAIN STEP****: START', time.asctime())

    tso_bamFolder = 'tso_bamFolder'
    import os
    tso_bamFolderPath = os.path.join(os.getcwd(), tso_bamFolder)
    os.mkdir(tso_bamFolderPath)
    print('tso_bamFolderPath', tso_bamFolderPath)

    outpath = args.out  # NEED A FRIENDLY WAY TO CHECK IF PATH EXISTS
    BAM = args.BAMFILE
    full_path_to_BAM = os.path.realpath(BAM)
    cores = args.cores

    # if outpath == None:
    #    outpath = os.path.dirname(full_path_to_BAM)

    # if not os.path.isdir(outpath):
    #    print('Warning: output path not found')

    #OUTBAM_FILENAME = os.path.join(outpath, 'TSO_reads.bam')
    #TMPBAM_FILENAME = os.path.join(outpath, 'tmp_tso.bam')

    #write_degraded_bam(full_path_to_BAM, TMPBAM_FILENAME, OUTBAM_FILENAME)

    #jsonFolder = 'jsonFolder'
    #jsonPath = os.path.join(os.getcwd(), jsonFolder)
    # os.mkdir(jsonPath)
    #print('jsonPath: ', jsonPath)

    print("Creating counter")
    counter = TSOReadCounter(args.BAMFILE, cores=args.cores, action="CountReads")
    print('****MAIN STEP****: FINISHED THE READCOUNTER', time.asctime())

    print("Registering files")
    counter.register_files()
    print('****MAIN STEP****: FINISHED REGISTERING', time.asctime())

    print("Counting reads: {}".format(args.BAMFILE))
    good_chromosomes = counter.run()
    print('****MAIN STEP****: FINISHED counter.run()', time.asctime())

    #print("Collecting read counts: {}".format(args.outdictionary))
    counter.combine(good_chromosomes, chrom_sep="\n")
    print('****MAIN STEP****: FINISHED COMBINING', time.asctime())

    print('Combine bam files', time.asctime())

    print("Merging bam files", time.asctime())
    bamFiles = glob.glob('tso_bamFolder/*.bam')
    print('bamFiles: ', bamFiles)
    #pysam.merge('-f', 'finalBam.bam', *bamFiles)
    pysam.cat(*bamFiles)
    print('****MAIN STEP****: FINISHED MERGING', time.asctime())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('--out', dest='out', help='Output folder', default=None)
    parser.add_argument('-C', '--cores', default=10, help="Number of cores.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    main(args)

# if __name__ == "__main__":
#    main(sys.argv[1:])
