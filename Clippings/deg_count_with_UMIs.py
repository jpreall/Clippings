#!/usr/bin/env python
# $ -cwd
# $ -v PATH,LD_LIBRARY_PATH
"""
Count mapped reads harboring TSO sub-sequene
Input:
    BAM file contains trimmed reads
    Only works on Cell Ranger v4 or newer output without pre-trimming TSO
Author: jpreall@cshl.edu
Date: 2021-04
"""
import argparse
import sys
import pysam
import collections
import sys
import os
import csv
import string
import argparse
import gzip
import scipy.sparse as sp
from scipy import io
import shutil
import numpy as np
import h5sparse
import h5py
import gtfparse
import time
from pararead import ParaReadProcessor
import logmuse
import json
import glob
import os.path as path
import pkg_resources
import collections
from collections import defaultdict


def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Degradation dictionary using ParaReadProcessor "
                    "implementation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "readsfile", help="Path to sequencing reads file (bam file).")

    parser.add_argument(
        "-D", "--outdictionary", help="Path to output dictionary.", default="dictionary.txt")

    parser.add_argument(
        "-C", "--cores", default=10, help="Number of cores.")

    parser.add_argument(
        "-TSS", "--TSSgtf", help="Path to gtf file.")

    parser.add_argument('--outdir', dest='outdir', help='Output folder', default="deg_raw_clipped")
    parser.add_argument('--genome', dest='genome',
                        help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
    parser.add_argument('--mtx', dest='mtx', help='Write output in 10X mtx format', default=False)
    parser.add_argument('--write_degraded_bam_file', dest='write_degraded_bam_file',
                        help='Writes all TSS filtered reads to a file called degraded_not_TSS.bam in parent directory', default=False)
    parser.add_argument('--include_introns', dest='include_introns',
                        action='store_true', help='Include intronic reads for nuclei data')

    parser = logmuse.add_logging_options(parser)
    return parser.parse_args(cmdl)


class ReadCounter(ParaReadProcessor):
    """Create child class ReadCounter from parent ParaReadProcessor.

    ReadCounter contains additional arguements such as feature dictionary
    and transcription start site dictionary.

    """

    def __init__(self, *args, **kwargs):
        """
        Extension of ParaReadProcessor class.

        Args:
            *args (list): ParaReadProcessor arguments.
            **kwargs (list): List containing additional 'features', '
            TSSdictionary', 'write_degraded_bam_file', and 'include_introns'
            arguments.

        """

        if 'features' in kwargs:
            features = kwargs.pop('features')
        else:
            args = list(features)
            features = args.pop()
        if 'TSSdictionary' in kwargs:
            TSSdictionary = kwargs.pop('TSSdictionary')
        else:
            TSSdictionary = args.pop()
        if 'write_degraded_bam_file' in kwargs:
            write_degraded_bam_file = kwargs.pop('write_degraded_bam_file')
        else:
            args = list(write_degraded_bam_file)
            write_degraded_bam_file = args.pop()
        if 'include_introns' in kwargs:
            include_introns = kwargs.pop('include_introns')
        else:
            args = list(include_introns)
            include_introns = args.pop()

        ParaReadProcessor.__init__(self, *args, **kwargs)
        self._features = features
        self._TSSdictionary = TSSdictionary
        self._write_degraded_bam_file = write_degraded_bam_file
        self._include_introns = include_introns

    def __call__(self, chromosome, _=None):
        """For each chromosome, perform a specific action such as building the
        degradation dictionary and writing out as a json file to be read in later.

        """

        # PYSAM CAN ONLY READ THROUGH A FILE ONCE. THEN YOU NEED TO RELOAD IT
        # THEREFORE IF DOING A LOOP USING PYSAM FUNCTIONS, DOING IT A SECOND
        # TIME WILL GIVE AN EMPTY OUTPUT

        reads = self.fetch_chunk(chromosome)

        deg_count_dict = bam_parser_chunk(reads, self._TSSdictionary, self._features,
                                          self._write_degraded_bam_file, self._include_introns)
        print('Finished running bam_parser_chunk: ', time.asctime())

        with open(self._tempf(chromosome), 'w') as f:
            #f.write("{}\t{}".format(chromosome, n_reads))
            print(chromosome, file=f)
            # the following might be too long of a print
            print(dict(list(deg_count_dict.items())[0:10]), file=f)

        print('Time Started for json write-out:', time.asctime())
        # 2021.08.09 Writing out dictionary as json to then re-load and combine
        jsonFolder = 'jsonFolder'
        import os
        jsonPath = os.path.join(os.getcwd(), jsonFolder)
        # change back if necessary
        filePath = os.path.join(jsonPath, str(chromosome)+'.json')
        with open(filePath, 'w') as f:
            json.dump(deg_count_dict, f)

        print('Finished json write-out:', time.asctime())
        return chromosome


def bam_parser_chunk(chunk, TSS_dict, feature_dictionary, write_degraded_bam_file, include_introns):
    """
    Go through every read and create a degradation dictionary of barcode and counts.

    For uniquely mapped, exonic reads with tso tags, get cell barcode, gene id and
    gene name tags. Creates a dictionary if cell barcode and degraded counts.

    Args:
        chunk (type): A chunk (chromosome) of sequencing reads from bam file.
        TSS_dict (dict): Dictionary of transcription start sites.
        feature_dictionary (dict): Dictionary of gene ids and gene names.
        write_degraded_bam_file (bool): True to write a degraded bam file. False to not.
        intron_introns (bool):

    Returns:
        dict: Dictionary of cell barcodes and degraded counts.

    """
    def tag_reader(aln, deg_count_dict, feature_dictionary):
        total_counts = 0
        total_TSS_counts = 0
        NO_UB = 0
        # Get Cell barcode and update set
        if aln.has_tag("CB"):
            cell_barcode = aln.get_tag("CB")
            if cell_barcode not in cell_barcode_set:
                cell_barcode_set.add(cell_barcode)
                deg_count_dict[cell_barcode] = collections.defaultdict(list)

            # Get Gene Name and Ensembl ID, if it there is one
            if aln.has_tag("GX"):
                gene_id = aln.get_tag("GX")

                # only take reads that unambiguously map to one gene
                VALID_GENE = gene_id in feature_dictionary.keys()

                if VALID_GENE:
                    TSSes_for_gene = TSS_dict[gene_id]
                    all_TSS_overlaps = set()

                    for TSS in TSSes_for_gene:
                        all_TSS_overlaps.add(aln.get_overlap(TSS-20, TSS+20))

                    NOT_TSS = max(all_TSS_overlaps) < 10
                    IS_TSS = not NOT_TSS

                    if aln.has_tag("GN"):
                        gene_name = aln.get_tag("GN")

                    # keep a running tally of total gene mapping counts:
                    total_counts += 1

                    # Only take degraded RNAs
                    #CLIPPED = aln.has_tag("ts:i")

                    #if NOT_TSS & VALID_GENE & CLIPPED:
                    #if NOT_TSS & VALID_GENE:
                    if VALID_GENE:
                        if gene_id not in deg_count_dict[cell_barcode].keys():
                            deg_count_dict[cell_barcode][gene_id] = set()

                        if aln.has_tag("UB"):
                            deg_count_dict[cell_barcode][gene_id].add(aln.get_tag("UB"))
                            if write_degraded_bam_file:
                                OUTPUT_BAM.write(aln)
                        else:
                            NO_UB += 1

                    # keep a running tally of all TSS mapping counts:
                   # elif IS_TSS & CLIPPED & VALID_GENE:
                    elif IS_TSS & VALID_GENE:
                        total_TSS_counts += 1

        return deg_count_dict, total_counts, total_TSS_counts, NO_UB

    # Initialize dictionary to store gene-cell matrix
    deg_count_dict = collections.defaultdict(list)

    # Maintain all detected barcodes as a set.  In future, this might be better to save the 10X whitelist in the repo and read from that.
    cell_barcode_set = set()

    # Running tallies of gene counts for QC, debugging, and reporting purposes
    total_lines_read = 0
    exon_count = 0
    intron_count = 0
    intergenic_count = 0
    exon_intron_count = 0

    if type(TSS_dict) != dict:
        print('Warning, TSS dictionary generation failed.  Exiting.')
        exit(0)

    #alignments = pysam.AlignmentFile(bamfile, "rb")

    # need to check to see if this feature works
    if write_degraded_bam_file:
        OUTPUT_BAM = pysam.AlignmentFile('degraded_not_TSS.bam', "wb", template=alignments)
        print('Writing TSS-filtered degraded BAM file to degraded_not_TSS.bam')

    print('Counting degraded reads...', time.asctime())

    print("include_introns is True", include_introns is True)
    if include_introns:
        for aln in chunk:
            #print("aln reached")
            total_lines_read += 1
            if total_lines_read % 1e7 == 0:
                print('Lines read:', f'{total_lines_read:,}')
            # Only consider uniquely mapped reads:
            if aln.mapping_quality == 255:
                #print("Mapping quality == 255")
                # Get both exonic and intronic reads
                if ((aln.get_tag("RE") == 'E') or (aln.get_tag("RE") == 'N')):
                    #print("Exon or intron")
                    if aln.get_tag("RE") == 'E':
                        exon_count += 1
                        # print("Exon")
                    elif aln.get_tag("RE") == 'N':
                        # print("Intron")
                        intron_count += 1
                    exon_intron_count += 1
                    #print("Tag reader start")
                    deg_count_dict, total_counts, total_TSS_counts, NO_UB = tag_reader(
                        aln, deg_count_dict, feature_dictionary)
                    #print("Tag reader end")
                elif aln.get_tag("RE") == 'I':
                    intergenic_count += 1
                    # print("Intergenic")
    else:
        print("include_introns is False", include_introns is False)
        for aln in chunk:
            #print("False version aln reached")
            total_lines_read += 1
            if total_lines_read % 1e7 == 0:
                print('Lines read:', f'{total_lines_read:,}')
            # Only consider uniquely mapped reads:
            #if aln.mapping_quality == 255:
            if aln.has_tag("xf:i"):
                #if (aln.get_tag("xf:i") & 8): # working, same as original
                if aln.get_tag("xf:i") == 17 or aln.get_tag("xf:i") == 25: #working, same as 8 after filter min_cells=1
                    #print(aln.get_tag("xf:i"))
                    # Only get exonic reads
                    if aln.get_tag("RE") == 'E':
                        #print("Exon only")
                        exon_count += 1
                        exon_intron_count += 1
                        deg_count_dict, total_counts, total_TSS_counts, NO_UB = tag_reader(
                            aln, deg_count_dict, feature_dictionary)
                    elif aln.get_tag("RE") == 'N':
                        # print("Intron")
                        intron_count += 1
                        exon_intron_count += 1
                    elif aln.get_tag("RE") == 'I':
                        # print("Intergenic")
                        intergenic_count += 1
    print('Total exonic: ', exon_count)
    print('Total intronic: ', intron_count)
    print('Total exon_intron: ', exon_intron_count)
    print('Total intergenic: ', intergenic_count)
    # broke for total_counts 2021.09.27 maybe because some iterations through aln do not have GN tag for gene name, thus creating an error(?)
    #print('Total counts:', total_counts)
    # broke: UnboundLocalError: local variable 'total_TSS_counts' referenced before assignment
    #print('Total TSS reads (not UMI de-duplicated):', total_TSS_counts)
    #print("Number of lines without 'UB' tag ", NO_UB)
    time.asctime()
    # Finally, collapse the UMI-containing sets into digital gene expression counts:
    print("Collapsing deg_count_dict to have the number of UMI", time.asctime())
    for cell in deg_count_dict.keys():
        for gene in deg_count_dict[cell].keys():
            deg_count_dict[cell][gene] = len(deg_count_dict[cell][gene])

    if write_degraded_bam_file:
        OUTPUT_BAM.close()

    return deg_count_dict  # , feature_dictionary


def merge_dicts(*dicts):
    """
    Merge each chromosome/chunk dictionary of barcodes, genes, and umi counts into a
    single degradation dictionary.

    Args:
        *dicts (tuple): Tuple of degradation dictionaries for each chunk/chromosome to be merged.

    Returns:
        dict: A merged dictionary containing all the keys and values for given tuple of dictionaries.

    """

    merged = defaultdict(dict)

    for d in dicts:
        for k, v in d.items():
            merged[k].update(v)

    return merged


def count_dict_to_sparse_matrix(data_dictionary, feature_dictionary):
    """
    write degraded counts to h5 or mtx format
    """
    barcodes = sorted(list(data_dictionary.keys()))
    bcdict = dict(zip(barcodes, np.arange(len(barcodes))))

    features = list(feature_dictionary.keys())
    featdict = dict(zip(features, np.arange(len(features))))

    # create matrix
    import scipy.sparse as sp

    # Create the sparse matrix
    data, row, col = [], [], []
    SHAPET = (len(features), len(barcodes))
    SHAPE = (len(barcodes), len(features))

    print('Building Dictionary of counts for matrix writing...', time.asctime())
    for cell in data_dictionary.keys():
        bcindex = bcdict[cell]
        newcols = [featdict[gene] for gene in data_dictionary[cell].keys()]
        row.extend([bcindex] * len(newcols))
        col.extend(newcols)
        data.extend(list(data_dictionary[cell].values()))

    print('Converting to scipy sparse matrix...', time.asctime())
    MATRIX = sp.csr_matrix((data, (row, col)), shape=SHAPE, dtype='int32')

    return MATRIX, barcodes


def write_10_mtx(data_dictionary, feature_dictionary):
    """
    Check to make sure this feature works
    write degraded counts to mtx format
    """

    # Prepare the matrix and cell barcodes
    MATRIX, barcodes = count_dict_to_sparse_matrix(data_dictionary, feature_dictionary)

    # write barcodes.tsv.gz
    #barcodes_file = gzip.open(os.path.join(output_folder,'barcodes.tsv.gz'), 'wt')
    barcodes_file = gzip.open('barcodes.tsv.gz', 'wt')
    writer = csv.writer(barcodes_file, delimiter='\t')
    for line in barcodes:
        writer.writerow([line])
    barcodes_file.close()

    # write features.tsv.gz
    features = list(feature_dictionary.keys())
    #features_file = gzip.open(os.path.join(output_folder,'features.tsv.gz'), 'wt')
    features_file = gzip.open('features.tsv.gz', 'wt')
    writer = csv.writer(features_file, delimiter='\t')
    for line in feature_dictionary:
        writer.writerow([line, feature_dictionary[line], 'Gene Expression'])
    features_file.close()

    # Write the matrix
    #io.mmwrite(os.path.join(output_folder,"matrix.mtx"), MATRIX)
    # with open(os.path.join(output_folder,"matrix.mtx"),'rb') as mtx_in:
    #    with gzip.open(os.path.join(output_folder,"matrix.mtx.gz"),'wb') as mtx_gz:
    #        shutil.copyfileobj(mtx_in, mtx_gz)
    #os.remove(os.path.join(output_folder,"matrix.mtx") )

    io.mmwrite('matrix.mtx', MATRIX)
    with open('matrix.mtx', 'rb') as mtx_in:
        with gzip.open('matrix.mtx.gz', 'wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove('matrix.mtx')


def write_10x_h5(data_dictionary, feature_dictionary, LIBRARY_ID=None, CHEMISTRY=None, genome=None):
    """
    write degraded counts to h5 format
    """
    # Prepare the matrix and cell barcodes
    MATRIX, barcodes = count_dict_to_sparse_matrix(data_dictionary, feature_dictionary)
    SHAPET = (len(feature_dictionary.keys()), len(barcodes))
    SHAPE = (len(barcodes), len(feature_dictionary.keys()))

    # Declare the output h5 file:
    outfile = 'raw_clipped_features_matrix.h5'
    print('Writing to '+outfile)

    # Encode Cell Barcodes
    BCS = np.array(barcodes).astype('S')

    # Encode Features
    FEATURES = np.array(list(feature_dictionary.values())).astype('S')
    FEATURE_IDS = np.array(list(feature_dictionary.keys())).astype('S')

    # May want to write a genome detection module in the future
    if genome:
        GENOME = genome
    else:
        print('No genome specified, writing attribute as unspecified_genome')
        GENOME = 'unspecified_genome'

    # Chemistry
    if not CHEMISTRY:
        print('No chemistry version specified, writing attribute as unspecified_chemistry')
        CHEMISTRY = 'unspecified_chemistry'

    # Sample name
    if not LIBRARY_ID:
        print('No library ID specified, writing attribute as unknown_library')
        LIBRARY_ID = 'unknown_library'

    # Other fields needed by Cellranger h5
    all_tag_keys = np.array([b'genome'])
    ORIG_GEM_GROUPS = np.array([1])

    # Write the h5
    print('Starting to write h5:', time.asctime())
    with h5sparse.File(outfile, 'w') as h5f:
        h5f.create_dataset('matrix/', data=MATRIX, compression="gzip")
        h5f.close()

    with h5py.File(outfile, 'r+') as f:
        f.create_dataset('matrix/barcodes', data=BCS)
        f.create_dataset('matrix/shape', (2,), dtype='int32', data=SHAPET)
        features = f.create_group('matrix/features')
        features.create_dataset('_all_tag_keys', (1,), 'S6', data=all_tag_keys)
        features.create_dataset('feature_type', data=np.array([b'Gene Expression'] * SHAPET[0]))
        features.create_dataset('genome', data=np.array([GENOME.encode()] * SHAPET[0]))
        features.create_dataset('id', data=FEATURE_IDS)
        features.create_dataset('name', data=FEATURES)

        f.attrs['chemistry_description'] = CHEMISTRY.encode()
        f.attrs['filetype'] = 'matrix'
        f.attrs['library_ids'] = LIBRARY_ID.encode()
        f.attrs['original_gem_groups'] = ORIG_GEM_GROUPS
        f.attrs['version'] = 2
        f.close()


def get_metadata(bamfile):
    """
    Read Library ID from bam file.
    """
    print("get pkg_resources: ", pkg_resources.resource_filename(
        'Clippings', 'files/737K-august-2016.txt.gz'))
    alignments = pysam.AlignmentFile(bamfile, "rb")

    # Need a backup if this step fails...
    # Is this even compatible with other versions of Cellranger?
    LIBRARY_ID = alignments.header.to_dict()['RG'][0]['SM']

    UMILENS = {}
    # Read the first 10000 lines and get UMI lengths
    for aln in alignments.head(10000):
        if aln.has_tag("CB"):
            if aln.has_tag("UB"):
                UMILEN = len(aln.get_tag("UB"))
                UMILENS[UMILEN] = UMILENS.get(UMILEN, 0) + 1

    PREDOMINANT_UMILEN = max(UMILENS, key=UMILENS.get)

    if PREDOMINANT_UMILEN == 12:
        CHEMISTRY = 'Single Cell 3\' v3'
    elif PREDOMINANT_UMILEN == 10:
        CHEMISTRY = 'Single Cell 3\' v2'
    else:
        CHEMISTRY = 'unspecified_chemistry'

    # Get the barcode whitelist for the relevenat chemistry
    BC_WHITELIST = fetch_barcode_whitelist(CHEMISTRY)

    return CHEMISTRY, LIBRARY_ID, BC_WHITELIST


def fetch_barcode_whitelist(CHEMISTRY):
    import pkg_resources
    VALID_CHEMISTRIES = {'Single Cell 3\' v2': pkg_resources.resource_filename('Clippings', 'files/737K-august-2016.txt.gz'),
                         'Single Cell 3\' v3': pkg_resources.resource_filename('Clippings', 'files/3M-february-2018.txt.gz'),
                         'unspecified_chemistry': pkg_resources.resource_filename('Clippings', 'files/3M-february-2018.txt.gz')
                         }

    WHITELIST_FILE = VALID_CHEMISTRIES[CHEMISTRY]

    if CHEMISTRY == 'unspecified_chemistry':
        print('Warning, unknown chemistry detected.  Defaulting to 10X Genomics Single Cell 3\' v3')

    BC_WHITELIST = set()
    with gzip.open(WHITELIST_FILE, mode='rb') as f:
        for line in f:
            BC_WHITELIST.add(line.decode('utf-8').strip())
    return BC_WHITELIST


def dict_of_TSSes(TSSgtf):

    print('Reading GTF file:', TSSgtf)
    report_time()
    df = gtfparse.read_gtf(TSSgtf)

    print('Building dictionary of valid gene_id:gene_name pairs...')
    report_time()
    feature_dictionary = df.loc[:, ['gene_id', 'gene_name']
                                ].drop_duplicates().set_index('gene_id').to_dict()['gene_name']

    print('Parsing out lines annotated as \'transcript\'...')
    report_time()
    transcript_df = df[df['feature'] == 'transcript'].loc[:, ['gene_name', 'gene_id',
                                                              'seqname', 'start', 'end', 'strand']].set_index('gene_id').drop_duplicates()

    print('Building list of valid TSSes per gene...')
    report_time()
    result = {}
    for gene_id, anno in transcript_df.iterrows():
        result[gene_id] = result.get(gene_id, set())
        result[gene_id].add(anno['start'] if anno['strand'] == '+' else anno['end'])
    return result, feature_dictionary


def report_time():
    print('Time started:', time.asctime())


def main(cmdl):
    """
    runner
    """
    # FileNotFoundError: [Errno 2] No such file or directory: '/cm/local/apps/uge/var/spool.p7444/bam08/files/3M-february-2018.txt.gz'
    FILES_PATH = path.abspath(path.join(path.dirname(__file__), "../files/"))
    print("This is the absolute file path: ", FILES_PATH)
    print("This is pkg_resources: ", pkg_resources.resource_filename(
        'Clippings', 'files/737K-august-2016.txt.gz'))

    import os
    args = _parse_cmdl(cmdl)
    global _LOGGER
    _LOGGER = logmuse.logger_via_cli(args, make_root=True)

    outdir = args.outdir
    genome = args.genome
    write_degraded_bam = args.write_degraded_bam_file

    assert not os.path.isdir(outdir), "Output directory already exists"

    os.mkdir(outdir)
    os.chdir(outdir)

    _LOGGER.debug("Run dict of TSSes")
    if args.TSSgtf != None:
        TSS_dict, feature_dictionary = dict_of_TSSes(args.TSSgtf)
    else:
        print("No TSS file supplied")

    # 2021.08.09 Creating directory to write out json dict
    jsonFolder = 'jsonFolder'
    jsonPath = os.path.join(os.getcwd(), jsonFolder)
    os.mkdir(jsonPath)
    print('jsonPath: ', jsonPath)

    _LOGGER.debug("Creating counter")
    counter = ReadCounter(args.readsfile, cores=args.cores,
                          outfile=args.outdictionary, action="CountReads",
                          TSSdictionary=TSS_dict, features=feature_dictionary,
                          write_degraded_bam_file=write_degraded_bam, include_introns=args.include_introns)
    print('****MAIN STEP****: FINISHED THE READCOUNTER', time.asctime())

    _LOGGER.debug("Registering files")
    counter.register_files()
    print('****MAIN STEP****: FINISHED REGISTERING', time.asctime())

    _LOGGER.info("Counting reads: {}".format(args.readsfile))
    good_chromosomes = counter.run()
    print('****MAIN STEP****: FINISHED counter.run()', time.asctime())

    _LOGGER.info("Collecting read counts: {}".format(args.outdictionary))
    counter.combine(good_chromosomes, chrom_sep="\n")
    print('****MAIN STEP****: FINISHED COMBINING', time.asctime())

    print('Time ended:', time.asctime())
    jsonFiles = glob.glob('./jsonFolder/*.json')
    print('jsonFiles: ', jsonFiles)

    # load dictionaries given name
    print('Load json Dicts time started:', time.asctime())
    jsonDictsList = []
    for json_entry in jsonFiles:
        with open(json_entry) as file:
            jsonDict = json.load(file)
        jsonDictsList.append(jsonDict)
    print('Finish load json Dicts time ended:', time.asctime())

    print('Merge dictionaries time started:', time.asctime())
    mergedDict = merge_dicts(*jsonDictsList)
    # print(dict(list(mergedDict.items())[0:15]))
    print('Merge dictionaries time ended:', time.asctime())

    print('Gathering metadata from bam file...')
    print('Time started:', time.asctime())
    CHEMISTRY, LIBRARY_ID, BC_WHITELIST = get_metadata(args.readsfile)

    # write 10X mtx format
    if args.mtx:
        try:
            print('Writing 10X-formatted mtx directory...', time.asctime())
            write_10_mtx(mergedDict, feature_dictionary)

        except IOError:
            print("I/O error")

    print('Writing 10X-formatted h5 file...', time.asctime())
    write_10x_h5(mergedDict, feature_dictionary, LIBRARY_ID, CHEMISTRY, genome=genome)

    print('Done!', time.asctime())


if __name__ == "__main__":
    # Override sys.argv
    #sys.argv = ['deg_count_with_UMIs.py', '/mnt/grid/scc/data/Preall/Preall_CR01/count/Preall_CR01_H_neg/outs/sorted_subsampledBAM.bam',
    #            '--TSSgtf', '/mnt/grid/scc/data/CellRanger/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf',
    #            '--outdir', 'testing_feb10', '--mtx', 'True']
    main(sys.argv[1:])
