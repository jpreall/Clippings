"""
Count mapped reads harboring TSO sub-sequene
Input:
    BAM file contains trimmed reads
    Only works on Cell Ranger v4 or newer output without pre-trimming TSO
Author: jpreall@cshl.edu
Date: 2021-04
"""
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

def bam_parser(bamfile, TSS_dict, feature_dictionary):

    # Initialize dictionary to store gene-cell matrix
    deg_count_dict = collections.defaultdict(list)

    # Maintain all detected barcodes as a set.  In future, this might be better to save the 10X whitelist in the repo and read from that.
    cell_barcode_set = set()

    # Running tallies of gene counts
    total_lines_read = 0
    total_counts = 0
    total_degraded_counts = 0
    total_TSS_counts = 0

    if type(TSS_dict) != dict:
        print('Warning, TSS dictionary generation failed.  Exiting.')
        exit(0)

    alignments = pysam.AlignmentFile(bamfile, "rb")
    print('Counting degraded reads...')
    for aln in alignments.fetch(until_eof=True):
    
        total_lines_read += 1 

        if total_lines_read % 1e7 == 0:
            print('Lines read:',f'{total_lines_read:,}')

        ### Only consider uniquely mapped reads:
        if aln.mapping_quality == 255:
        
            ### Get Cell barcode and update set
            if aln.has_tag("CB"):
                cell_barcode = aln.get_tag("CB")
                if cell_barcode not in cell_barcode_set:
                    cell_barcode_set.add(cell_barcode)
                    deg_count_dict[cell_barcode] = collections.defaultdict(list)

                ### Get Gene Name and Ensembl ID, if it there is one
                if aln.has_tag("GX"):
                    gene_id = aln.get_tag("GX")

                    ### only take reads that unambiguously map to one gene
                    VALID_GENE = gene_id in feature_dictionary.keys()

                    if VALID_GENE:
                        TSSes_for_gene = TSS_dict[gene_id]
                        all_TSS_overlaps = set()

                        for TSS in TSSes_for_gene:
                            all_TSS_overlaps.add(aln.get_overlap(TSS-20,TSS+20))

                        NOT_TSS = max(all_TSS_overlaps) < 10
                        IS_TSS = not NOT_TSS
                        

                        if aln.has_tag("GN"):
                            gene_name = aln.get_tag("GN")

                        #### keep a running tally of total gene mapping counts:
                        total_counts += 1

                        #### Only take degraded RNAs
                        CLIPPED = aln.has_tag("ts:i")

                        if NOT_TSS & VALID_GENE & CLIPPED:
                            if gene_id not in deg_count_dict[cell_barcode].keys():
                                deg_count_dict[cell_barcode][gene_id] = set()

                            NO_UB = 0
                            if aln.has_tag("UB"):
                                deg_count_dict[cell_barcode][gene_id].add(aln.get_tag("UB"))
                            else:
                                NO_UB += 1
                                
                        #keep a running tally of all TSS mapping counts:
                        elif IS_TSS & CLIPPED & VALID_GENE:
                            total_TSS_counts += 1

            else:
                continue

    #print('Total counts:', total_counts)
    #print('Total degraded counts:', total_degraded_counts)
    print('Total TSS reads (not UMI de-duplicated):', total_TSS_counts)
    print("Number of lines without 'UB' tag ", NO_UB)
        
    ## Finally, collapse the UMI-containing sets into digital gene expression counts:		
    for cell in deg_count_dict.keys():
        for gene in deg_count_dict[cell].keys():
            deg_count_dict[cell][gene] = len(deg_count_dict[cell][gene])

    return deg_count_dict, feature_dictionary
    

def count_dict_to_sparse_matrix(data_dictionary, feature_dictionary):
    """
    write degraded counts to h5 or mtx format
    """
    barcodes = sorted(list(data_dictionary.keys()))
    bcdict = dict(zip(barcodes,np.arange(len(barcodes))))

    features = list(feature_dictionary.keys())
    featdict = dict(zip(features,np.arange(len(features))))

    ## create matrix
    import scipy.sparse as sp

    # Create the sparse matrix
    data, row, col = [], [], []
    SHAPET = (len(features),len(barcodes))
    SHAPE = (len(barcodes),len(features))

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
    
    
def write_10_mtx(data_dictionary, feature_dictionary, output_folder):
    """
    write degraded counts to mtx format
    """
    
    #Prepare the matrix and cell barcodes
    MATRIX, barcodes = count_dict_to_sparse_matrix(data_dictionary, feature_dictionary)

    ## write barcodes.tsv.gz
    barcodes_file = gzip.open(os.path.join(output_folder,'barcodes.tsv.gz'), 'wt')
    writer = csv.writer(barcodes_file, delimiter='\t')
    for line in barcodes:
        writer.writerow([line])
    barcodes_file.close()

    ## write features.tsv.gz
    features = list(feature_dictionary.keys())
    features_file = gzip.open(os.path.join(output_folder,'features.tsv.gz'), 'wt')
    writer = csv.writer(features_file, delimiter='\t')
    for line in feature_dictionary:
        writer.writerow([line,feature_dictionary[line],'Gene Expression'])
    features_file.close()

 	# Write the matrix
    io.mmwrite(os.path.join(output_folder,"matrix.mtx"), MATRIX)
    with open(os.path.join(output_folder,"matrix.mtx"),'rb') as mtx_in:
        with gzip.open(os.path.join(output_folder,"matrix.mtx.gz"),'wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(os.path.join(output_folder,"matrix.mtx") )

def write_10x_h5(data_dictionary, feature_dictionary, output_folder, LIBRARY_ID=None, CHEMISTRY=None, genome=None):
    """
    write degraded counts to h5 format
    """
    #Prepare the matrix and cell barcodes
    MATRIX, barcodes = count_dict_to_sparse_matrix(data_dictionary, feature_dictionary)
    SHAPET = (len(feature_dictionary.keys()),len(barcodes))
    SHAPE = (len(barcodes),len(feature_dictionary.keys()))
    
    #Declare the output h5 file:
    outfile=os.path.join(output_folder,'raw_clipped_features_matrix.h5')
    print('Writing to '+outfile)

    # Encode Cell Barcodes
    BCS = np.array(barcodes).astype('S')

    # Encode Features
    FEATURES = np.array(list(feature_dictionary.values())).astype('S')
    FEATURE_IDS = np.array(list(feature_dictionary.keys())).astype('S')

    ## May want to write a genome detection module in the future
    if genome:
        GENOME=genome
    else:
        print('No genome specified, writing attribute as unspecified_genome')
        GENOME='unspecified_genome'

    ## Chemistry
    if not CHEMISTRY:
        print('No chemistry version specified, writing attribute as unspecified_chemistry')
        CHEMISTRY='unspecified_chemistry'

    # Sample name
    if not LIBRARY_ID:
        print('No library ID specified, writing attribute as unknown_library')
        LIBRARY_ID = 'unknown_library'

    # Other fields needed by Cellranger h5
    all_tag_keys = np.array([b'genome'])
    ORIG_GEM_GROUPS = np.array([1])

    # Write the h5
    print('Starting to write h5:',time.asctime())
    with h5sparse.File(outfile, 'w') as h5f:
        h5f.create_dataset('matrix/', data=MATRIX, compression="gzip")
        h5f.close()

    with h5py.File(outfile, 'r+') as f:
        f.create_dataset('matrix/barcodes', data=BCS)
        f.create_dataset('matrix/shape', (2,),dtype='int32', data=SHAPET)
        features = f.create_group('matrix/features')
        features.create_dataset('_all_tag_keys', (1,),'S6', data=all_tag_keys)
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
    alignments = pysam.AlignmentFile(bamfile, "rb")

    # Need a backup if this step fails...
    # Is this even compatible with other versions of Cellranger?
    LIBRARY_ID = alignments.header.to_dict()['RG'][0]['SM']

    UMILENS = {}
    # Read the first 10000 lines and get UMI lengths
    for aln in alignments.head(10000):
                if aln.has_tag("CB"):
                    if aln.has_tag("UB"):
                        UMILEN=len(aln.get_tag("UB"))
                        UMILENS[UMILEN] = UMILENS.get(UMILEN, 0) + 1

    PREDOMINANT_UMILEN = max(UMILENS, key=UMILENS.get)

    if PREDOMINANT_UMILEN == 12:
        CHEMISTRY = 'Single Cell 3\' v3'
    elif PREDOMINANT_UMILEN == 10:
        CHEMISTRY = 'Single Cell 3\' v2'
    else:
        CHEMISTRY = 'unspecified_chemistry'

    ## Get the barcode whitelist for the relevenat chemistry
    BC_WHITELIST = fetch_barcode_whitelist(CHEMISTRY)

    return CHEMISTRY, LIBRARY_ID, BC_WHITELIST



def fetch_barcode_whitelist(CHEMISTRY):
    VALID_CHEMISTRIES = {
        'Single Cell 3\' v2':'../files/737K-august-2016.txt.gz',
        'Single Cell 3\' v3':'../files/3M-february-2018.txt.gz',
        'unspecified_chemistry':'../files/3M-february-2018.txt.gz'
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
    feature_dictionary = df.loc[:,['gene_id','gene_name']].drop_duplicates().set_index('gene_id').to_dict()['gene_name']

    print('Parsing out lines annotated as \'transcript\'...')
    report_time()
    transcript_df = df[df['feature'] == 'transcript'].loc[:,['gene_name','gene_id','seqname','start','end','strand']].set_index('gene_id').drop_duplicates()

    print('Building list of valid TSSes per gene...')
    report_time()
    result = {}
    for gene_id, anno in transcript_df.iterrows():
        result[gene_id] = result.get(gene_id, set())
        result[gene_id].add(anno['start'] if anno['strand'] == '+' else anno['end'])
    return result, feature_dictionary

def report_time():
    print('Time started:',time.asctime())
	
def main(args):
    outdir = args.out
    genome = args.genome

    if os.path.isdir(outdir):
        overwrite = input('\nOutput directory already exists. Overwrite? Y/N ')
        if overwrite.lower() == 'n':
            exit(0)
        elif overwrite.lower() == 'y':
            shutil.rmtree(outdir)

    os.mkdir(outdir)

    print('Gathering metadata from bam file...')
    print('Time started:',time.asctime())
    CHEMISTRY, LIBRARY_ID, BC_WHITELIST = get_metadata(args.bamfile)

    if args.TSSgtf != None:
        TSS_dict, feature_dictionary = dict_of_TSSes(args.TSSgtf)
        deg_count_dict, feature_dictionary = bam_parser(args.bamfile, TSS_dict, feature_dictionary)
    else:
        deg_count_dict, feature_dictionary = bam_parser_noTSS(args.bamfile)


    #write 10X mtx format
    if args.mtx:
        try:
            print('Writing 10X-formatted mtx directory...')
            write_10_mtx(deg_count_dict, feature_dictionary, outdir)

        except IOError:
            print("I/O error")

    ##write 10X h5 format
    try:
        print('Writing 10X-formatted h5 file...')
        write_10x_h5(deg_count_dict, feature_dictionary, outdir, LIBRARY_ID, CHEMISTRY, genome=genome)

    except IOError:
        print("I/O error")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile', help='Input bam file')
    parser.add_argument('--out', dest='out', help='Output folder', default="raw_clipped_features_matrix")
    parser.add_argument('--genome', dest='genome', help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
    parser.add_argument('--mtx', dest='mtx', help='Write output in 10X mtx format', default=False)
    parser.add_argument('--TSSgtf', dest='TSSgtf', help='GTF file for filtering TSS-overlapping reads', default=None)
    parser.add_argument('--miRNAgtf', dest='miRNAgtf', help='GTF file for assigning TSO-reads to miRNA cropping sites', default=None)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    main(args)
