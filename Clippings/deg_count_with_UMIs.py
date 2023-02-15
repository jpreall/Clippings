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
import gzip
import shutil
#import scipy.sparse as sp
from scipy import io
import numpy as np
import h5sparse
import h5py
import gtfparse
from multiprocessing import Pool
import logging
import json
import glob
import os.path as path
import pkg_resources

def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Degradation dictionary using ParaReadProcessor "
                    "implementation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "readsfile", help="Path to sequencing reads file (bam file).")

    #parser.add_argument(
    #    "-D", "--outdictionary", help="Path to output dictionary.", default="dictionary.txt")

    parser.add_argument(
        "-C", "--cores", default=10, help="Number of cores.")
    
    parser.add_argument(
        "-matrix_folder", "--matrix_folder", required=True, help='10x Genomics matrix folder (raw or filtered) to concatenate miRNAs into')

    parser.add_argument(
        "-TSS", "--TSSgtf", help="Path to gtf file.")

    parser.add_argument('--outdir', dest='outdir', help='Output folder', default="deg_raw_clipped")
    parser.add_argument('--genome', dest='genome',
                        help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
    parser.add_argument('--write_bam_output', dest='write_degraded_bam_file',
                        help='Writes all TSS filtered reads to a file called degraded_not_TSS.bam in parent directory', default=False)

    return parser.parse_args(cmdl)


def count_reads(dict_of_input):
    """
    Function to parallelize the downsampling by chrom
    
    Args:
        dict_of_input (dict): Dictionary containing all args being passed from parse bam file, includes "args", "TSS_dict", "feature_dictionary" and "jsonPath"
    
    """
    args = dict_of_input["args"]
    with pysam.AlignmentFile(args.readsfile, "rb") as bam:
        chroms = bam.references
    input_things = [(chrom, dict_of_input) for chrom in chroms]
    with Pool(int(args.cores)) as p:
        completed_chroms = p.starmap(bam_parser_chunk, input_things)

def parse_bam_file(args, TSS_dict, feature_dictionary):
    """
    Outer function to parallelize without pararead

    Args:
        args: passed in command line arguments
        TSS_dict (dict): Dictionary of transcription start sites.
        feature_dictionary (dict): Dictionary of gene ids and gene names.

    """
    #first make temp json folder 
    jsonFolder = 'jsonFolder'
    jsonPath = os.path.join(os.getcwd(), jsonFolder)
    count_reads({"args": args, "TSS_dict": TSS_dict, "feature_dictionary": feature_dictionary, "jsonPath": jsonPath})
    return jsonPath 


def bam_parser_chunk(chrom, dict_of_input):
    """
    Go through every read and create a degradation dictionary of barcode and counts.

    For uniquely mapped, reads with tso tags, get cell barcode, gene id and
    gene name tags. Creates a dictionary if cell barcode and degraded counts and writes it to json.

    Args:
        chrom (type): A chromosome from which we'll pull the reads from the bam file.
        dict_of_input (dict): Dictionary containing all args being passed from parse bam file, includes "args", "TSS_dict", "feature_dictionary" and "jsonPath".

    """
    args = dict_of_input['args']
    TSS_dict = dict_of_input['TSS_dict']
    feature_dictionary = dict_of_input['feature_dictionary']
    jsonPath = dict_of_input['jsonPath']

    write_degraded_bam_file = args.write_degraded_bam_file

    # Initialize dictionary to store gene-cell matrix
    deg_count_dict = collections.defaultdict(list)

    # Maintain all detected barcodes as a set.  In future, this might be better to save the 10X whitelist in the repo and read from that.
    cell_barcode_set = set()

    # Running tallies of gene counts for QC, debugging, and reporting purposes
    total_lines_read = 0
    total_counts = 0
    total_TSS_counts = 0

    if type(TSS_dict) != dict:
        logging.error('Warning, TSS dictionary generation failed.  Exiting.')
        # should probably be exit(1)
        exit(0)

    if write_degraded_bam_file:
        alignments = pysam.AlignmentFile(args.readsfile, "rb")
        OUTPUT_BAM = pysam.AlignmentFile(f"temp_bams/{chrom}_degraded_not_TSS.bam", "wb", template=alignments)
        #logging.info('Writing TSS-filtered degraded BAM file to degraded_not_TSS.bam')

    #logging.info(f'Counting degraded reads for {chrom}')
    bamin = pysam.AlignmentFile(args.readsfile, "rb")
    for aln in bamin.fetch(chrom):
        total_lines_read += 1
        if total_lines_read % 1e7 == 0:
            logging.info(f'Lines read: {total_lines_read:,} on {chrom}')
            
        #Skip any reads without all of the tags we need

        if not (aln.has_tag("CB")) & (aln.has_tag("GX")) & aln.has_tag("UB"):
            continue

        #Skip any reads that Cellranger didn't count as a UMI or possible UMI

        if not (aln.get_tag("xf") == 17) | (aln.get_tag("xf") == 25):
            continue

        # Get Cell barcode and update set
        
        cell_barcode = aln.get_tag("CB")
        if cell_barcode not in cell_barcode_set:
            cell_barcode_set.add(cell_barcode)
            deg_count_dict[cell_barcode] = collections.defaultdict(list)

        # Get Gene Name and Ensembl ID, if it there is one
        
        gene_id = aln.get_tag("GX")

        # only take reads that unambiguously map to one gene

        VALID_GENE = gene_id in feature_dictionary.keys()

        if not VALID_GENE:
            continue
    
        TSSes_for_gene = TSS_dict[gene_id]
        all_TSS_overlaps = set()

        for TSS in TSSes_for_gene:
            all_TSS_overlaps.add(aln.get_overlap(TSS-20, TSS+20))

        NOT_TSS = max(all_TSS_overlaps) < 10
        IS_TSS = not NOT_TSS

        # keep a running tally of total gene mapping counts:
        total_counts += 1

        # Only take degraded RNAs
        CLIPPED = aln.has_tag("ts:i")

        if NOT_TSS & CLIPPED:
            if gene_id not in deg_count_dict[cell_barcode].keys():
                deg_count_dict[cell_barcode][gene_id] = set()

            deg_count_dict[cell_barcode][gene_id].add(aln.get_tag("UB"))
            if write_degraded_bam_file:
                OUTPUT_BAM.write(aln)

        # keep a running tally of all TSS mapping counts:
        elif IS_TSS & CLIPPED:
            total_TSS_counts += 1

    # Finally, collapse the UMI-containing sets into digital gene expression counts:
    for cell in deg_count_dict.keys():
        for gene in deg_count_dict[cell].keys():
            deg_count_dict[cell][gene] = len(deg_count_dict[cell][gene])

    if write_degraded_bam_file:
        OUTPUT_BAM.close()

    filePath = os.path.join(jsonPath, str(chrom)+'.json')
    with open(filePath, 'w') as f:
        json.dump(deg_count_dict, f)
    logging.info(f"Finished writing deg_count_dict for {chrom}")

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

    merged = collections.defaultdict(dict)

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

    logging.info('Building Dictionary of counts for matrix writing...')
    for cell in data_dictionary.keys():
        bcindex = bcdict[cell]
        newcols = [featdict[gene] for gene in data_dictionary[cell].keys()]
        row.extend([bcindex] * len(newcols))
        col.extend(newcols)
        data.extend(list(data_dictionary[cell].values()))

    logging.info('Converting to scipy sparse matrix...')
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
    logging.info('Writing to '+outfile)

    # Encode Cell Barcodes
    BCS = np.array(barcodes).astype('S')

    # Encode Features
    FEATURES = np.array(list(feature_dictionary.values())).astype('S')
    FEATURE_IDS = np.array(list(feature_dictionary.keys())).astype('S')

    # May want to write a genome detection module in the future
    if genome:
        GENOME = genome
    else:
        logging.info('No genome specified, writing attribute as unspecified_genome')
        GENOME = 'unspecified_genome'

    # Chemistry
    if not CHEMISTRY:
        logging.info('No chemistry version specified, writing attribute as unspecified_chemistry')
        CHEMISTRY = 'unspecified_chemistry'

    # Sample name
    if not LIBRARY_ID:
        logging.info('No library ID specified, writing attribute as unknown_library')
        LIBRARY_ID = 'unknown_library'

    # Other fields needed by Cellranger h5
    all_tag_keys = np.array([b'genome'])
    ORIG_GEM_GROUPS = np.array([1])

    # Write the h5
    logging.info('Starting to write h5')
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
    logging.info(f"get pkg_resources: {pkg_resources.resource_filename('Clippings', 'files/737K-august-2016.txt.gz')}")
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
        logging.info('Warning, unknown chemistry detected.  Defaulting to 10X Genomics Single Cell 3\' v3')

    BC_WHITELIST = set()
    with gzip.open(WHITELIST_FILE, mode='rb') as f:
        for line in f:
            BC_WHITELIST.add(line.decode('utf-8').strip())
    return BC_WHITELIST


def dict_of_TSSes(TSSgtf):

    logging.info('Reading GTF file:')
    df = gtfparse.read_gtf(TSSgtf)

    logging.info('Building dictionary of valid gene_id:gene_name pairs...')
    feature_dictionary = df.loc[:, ['gene_id', 'gene_name']
                                ].drop_duplicates().set_index('gene_id').to_dict()['gene_name']

    logging.info('Parsing out lines annotated as \'transcript\'...')
    transcript_df = df[df['feature'] == 'transcript'].loc[:, ['gene_name', 'gene_id',
                                                              'seqname', 'start', 'end', 'strand']].set_index('gene_id').drop_duplicates()

    logging.info('Building list of valid TSSes per gene...')
    result = {}
    for gene_id, anno in transcript_df.iterrows():
        result[gene_id] = result.get(gene_id, set())
        result[gene_id].add(anno['start'] if anno['strand'] == '+' else anno['end'])
    return result, feature_dictionary


def main(cmdl):
    """
    runner
    """
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S",
                        level=logging.INFO, force=True)

    # FileNotFoundError: [Errno 2] No such file or directory: '/cm/local/apps/uge/var/spool.p7444/bam08/files/3M-february-2018.txt.gz'
    FILES_PATH = path.abspath(path.join(path.dirname(__file__), "../files/"))
    #logging.info("This is the absolute file path: "+FILES_PATH)
    #logging.info("This is pkg_resources: " + pkg_resources.resource_filename(
    #    'Clippings', 'files/737K-august-2016.txt.gz'))


    args = _parse_cmdl(cmdl)

    outdir = args.outdir
    genome = args.genome
    write_degraded_bam = args.write_degraded_bam_file

    assert not os.path.isdir(outdir), "Output directory already exists"

    os.mkdir(outdir)
    os.chdir(outdir)

    logging.debug("Run dict of TSSes")
    if args.TSSgtf != None:
        TSS_dict, feature_dictionary = dict_of_TSSes(args.TSSgtf)
        #deg_count_dict, feature_dictionary = bam_parser(args.bamfile, TSS_dict, feature_dictionary, write_degraded_bam_file)
    else:
        #deg_count_dict, feature_dictionary = bam_parser_noTSS(args.bamfile)
        logging.info("No TSS file supplied")

    # 2021.08.09 Creating directory to write out json dict
    jsonFolder = 'jsonFolder'
    jsonPath = os.path.join(os.getcwd(), jsonFolder)
    os.mkdir(jsonPath)
    logging.info('jsonPath: '+ jsonPath)
    if write_degraded_bam:
        bams_folder = os.path.join(os.getcwd(), "temp_bams")
        os.mkdir(bams_folder)
    
    #Parallel processing of the bam file to count degraded read
    parse_bam_file(args, TSS_dict, feature_dictionary)
    logging.info('****MAIN STEP****: FINISHED PARSING BAM')

    if write_degraded_bam:
        # steps to merge, sort, index and delete temporary bam files to write the degraded bam
        logging.info("Combining Degraded Bam Files")
        temp_bams = glob.glob(bams_folder+"/*.bam")
        print(temp_bams[0])
        pysam.merge(*["degraded_not_TSS_unsorted.bam"] + temp_bams)
        logging.info("Sorting Degraded Bam File")
        pysam.sort("-o", "degraded_not_TSS.bam", "degraded_not_TSS_unsorted.bam")
        os.remove("degraded_not_TSS_unsorted.bam")
        logging.info("Indexing Degraded Bam File")
        pysam.index("degraded_not_TSS.bam")
        [os.remove(i) for i in glob.glob(bams_folder+"/*.bam")]
        os.rmdir(bams_folder)
        logging.info("Degraded Bam Written")

    # Merge the jsons together into one dictionary 
    jsonFiles = glob.glob('./jsonFolder/*.json')

    logging.info('Load json Dicts')
    jsonDictsList = []
    for json_entry in jsonFiles:
        with open(json_entry) as file:
            jsonDict = json.load(file)
        jsonDictsList.append(jsonDict)
    logging.info('Finished load json Dicts')

    logging.info('Merge dictionaries')
    merged_unfiltered_Dict = merge_dicts(*jsonDictsList)

    #filtering to match supplied feature matrix, then cleaning up jsonFolder
    with gzip.open(args.matrix_folder+"/barcodes.tsv.gz") as f:
        raw_barcodes = f.readlines()
    barcode_list_from_matrix = [str(i).replace("\\n'", "").replace("b'", "") for i in raw_barcodes]
    mergedDict = {i:merged_unfiltered_Dict.get(i, collections.defaultdict(list)) for i in barcode_list_from_matrix}
    [os.remove(i) for i in jsonFiles]
    os.rmdir("jsonFolder")

    logging.info('Gathering metadata from bam file...')
    CHEMISTRY, LIBRARY_ID, BC_WHITELIST = get_metadata(args.readsfile)

    # write 10X mtx format
    try:
        logging.info('Writing 10X-formatted mtx directory...')
        write_10_mtx(mergedDict, feature_dictionary)

    except IOError:
        logging.error("I/O error")

    logging.info('Writing 10X-formatted h5 file...')
    write_10x_h5(mergedDict, feature_dictionary, LIBRARY_ID, CHEMISTRY, genome=genome)

    logging.info('Done!')


if __name__ == "__main__":
    main(sys.argv[1:])


