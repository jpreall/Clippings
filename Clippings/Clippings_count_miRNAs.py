#!/usr/bin/env python
#$ -cwd
#$ -v PATH,LD_LIBRARY_PATH
"""
Clippings_count_miRNAs.py

Count mapped reads harboring TSO sub-subsequence. Only works with CellRanger 4.0+ outputs.


Inputs:
    BAMFILE = BAM file produced by Cellranger v4+ with TSO reads tagged with ts:i
    REFERENCE_FILE = gff3 file with miRNA coordinates downloaded from miRBase
    ### TODO: Add ability to use miRNA GTF/coordinate files produced elsewhere

    Usage: Clippings_count_miRNA.py BAMFILE REFERENCE_FILE

Author: jpreall@cshl.edu
Date: 2021-07
"""

import sys
import numpy as np
import pandas as pd
import os
import glob
import pysam
import anndata
import argparse
import shutil
import h5sparse
import h5py
import gtfparse
import time
import scanpy as sc

##packages to try to do without:
#import scanpy as sc
#import scipy
#import seaborn as sns
#import matplotlib.pyplot as pl
#import gzip
#import csv
#import string

def attributes_to_columns_from_mirbase_GFF3(file):
    """
    Given a GFF3 file from mirbase, turn attributes into pandas dataframe columns.

    Args:
        file (string path): Path to .gff3 formatted file with miRNA annotations downloaded from miRBase.

    Returns:
        pandas DataFrame: Dataframe with columns containing the attributes from GFF3 file.

    """

    print("Starting attirubtes_to_columns_from_mirbase_GFF3: ", time.asctime())

    import pandas as pd
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    attribute_df = pd.read_table(file, comment='#', header=None)
    attribute_df.columns = columns

    #attribute_df = self.copy()
    #miRNA_anno_df = attribute_df.loc[:, "seqname":"attributes"]
    miRNA_anno_df = attribute_df.copy()

    attribute_df["at_dic"] = attribute_df.attributes.apply(
        lambda attributes: dict(
            [
                key_value_pair.split(sep="=", maxsplit=1)
                for key_value_pair in attributes.split(";")
            ]
        )
    )
    attribute_df["at_dic_keys"] = attribute_df["at_dic"].apply(
        lambda at_dic: list(at_dic.keys())
    )
    merged_attribute_list = []
    for key in attribute_df["at_dic_keys"]:
        merged_attribute_list += key

    nonredundant_list = sorted(list(set(merged_attribute_list)))
    for atr in nonredundant_list:
        miRNA_anno_df[atr] = attribute_df["at_dic"].apply(
            lambda at_dic: at_dic.get(atr)
        )
    print("Finish attirubtes_to_columns_from_mirbase_GFF3: ", time.asctime())
    return miRNA_anno_df

def dict_of_parent_names(miRNA_anno_df):
    """Extracts primary transcripts and returns a dictionary of ID and Name"""
    PARENT_DF = miRNA_anno_df[miRNA_anno_df['feature'] == 'miRNA_primary_transcript']
    PARENT_DICT = dict(zip(PARENT_DF['ID'],PARENT_DF['Name']))
    return PARENT_DICT


def make_dictionary_of_miRNA_Drosha_coords(miRNA_anno_df, PARENT_DICT):
    """
    Creates dictionary of chromosome to miRNA and it's drosha coordinates.

    Args:
        miRNA_anno_df (pandas DataFrame): Dataframe with columns containing the attributes from GFF3 file.
        PARENT_DICT (dict): Dictionary of ID to Name of only primary transcripts.

    Returns:
        dict: Dictionary of chromosome as key and value as default dict.
        Default dict has miRNA name and drosha coordinate with '+' or '-' strand.

    """

    #Only consider entries with a clearly marked 3p arm:
    print("Starting make_dictionary_of_miRNA_Drosha_coords: ", time.asctime())

    threep = miRNA_anno_df[miRNA_anno_df['Name'].str.match('.*3p$')]

    #Drop any entries where the 3p arm is erroneously marked a
    dupthreep = threep[threep['Derives_from'].duplicated()]
    threep = threep[~threep['Derives_from'].duplicated()]


    coord_dict = {}

    for __, line in threep.iterrows():
        COORD = line['start'] if line['strand'] == '-' else line['end']
        CHROM = line['seqname']
        DROSHA_SITE = CHROM+':' + str(COORD)
        DERIVES_FROM = line['Derives_from']
        coord_dict[DERIVES_FROM] = DROSHA_SITE

    coord_dict = {PARENT_DICT[DERIVES_FROM]:coord_dict[DERIVES_FROM] for DERIVES_FROM in coord_dict}

    import collections

    chrom_dict = {}
    for __, line in threep.iterrows():
        COORD = line['start'] if line['strand'] == '-' else line['end']
        CHROM = line['seqname']
        STRAND = line['strand']
        DROSHA_SITE = CHROM+':' + str(COORD)
        DERIVES_FROM = line['Derives_from']
        if CHROM not in chrom_dict.keys():
            chrom_dict[CHROM] = collections.defaultdict(list)
        chrom_dict[CHROM][PARENT_DICT[DERIVES_FROM]] = (COORD,STRAND)
    print("Finish make_dictionary_of_miRNA_Drosha_coords: ", time.asctime())
    return chrom_dict

def fix_chr_chromnames(chrom_dict, BAM):
    """
    Checks if BAM file is using 'chr' prefix on chromosome names and fixes the dictionary if necessary.

    Args:
        chrom_dict (dict): Dictionary of chromosome as key and value as default dict.
            Default dict has miRNA name and drosha coordinate with '+' or '-' strand.
        BAM (bamfile): Bamfile produced by 10x Genomics' CellRanger

    Returns:
        dict: Updated version of chrom_dict with correct 'chr' prefix.

    """

    print("Starting fix_chr_chromnames: ", time.asctime())
    bamfile = pysam.AlignmentFile(BAM, "rb")

    import re
    pattern = re.compile(r'chr[0-9]+', re.IGNORECASE)
    REF_STARTING_WITH_CHR = list(filter(pattern.match, bamfile.header.references))

    # only fix if < 10% of chromosome names don't start with chr
    if len(REF_STARTING_WITH_CHR) < 0.1*len(bamfile.header.references):
        newkeys = [name.replace('chr','') for name in chrom_dict.keys()]
        newkeys_pass_filter = [name for name in newkeys if name in bamfile.header.references]
        print('Changing',len(newkeys_pass_filter),'chromosome names in GTF file to match BAM file...')
        for newkey in newkeys_pass_filter:
            chrom_dict[newkey] = chrom_dict.pop('chr' + newkey)
    print("Finish fix_chr_chromnames: ", time.asctime())
    return chrom_dict



def count_miRNAs(BAM, chrom_dict, flanks=0):
    """
    Reads bamfile and stores read information into a pandas DataFrame and counts number
    of each miRNA.

    Args:
        BAM (bamfile): Bamfile produced by 10x Genomics' CellRanger
        chrom_dict (dict): Dictionary of chromosome as key and value as default dict.
            Default dict has miRNA name and drosha coordinate with '+' or '-' strand.
        flanks (int): Default 0.

    Returns:
        (tuple): tuple containing:

            count_table (pandas DataFrame): Dataframe of counts of each miRNA.
            example_read (pysam AlignmentSegment): Read as AlignmentSegment.
            results (pandas DataFrame): Dataframe of all reads and info from pysam.

    """

    print("Starting count_miRNAs: ", time.asctime())
    import collections
    example_read = None

    allmi = {}
    for c in chrom_dict.keys():
        for m in chrom_dict[c].keys():
            allmi[m] = c

    OUTBAM_FILENAME = 'miRNA_near_matches.bam'
    print(BAM)
    alignments = pysam.AlignmentFile(BAM, "rb")
    MI_READS = pysam.AlignmentFile('tmp_miRNA_matching.bam', "wb", template=alignments)


    count_table = collections.defaultdict(set)
    detailed_results = collections.defaultdict(list)
    tally = 0
    NO_UB = 0


    distances = []
    for i in range(len(allmi.keys())):
        m = list(allmi.keys())[i]
        chroms = [allmi[m]]
        mirnas = [m]

        for chrom in chroms:

            for mirna in mirnas:

                mirna_strand = chrom_dict[chrom][mirna][1]
                DROSHA_SITE = chrom_dict[chrom][mirna][0]

                if mirna_strand == '+':
                    start = chrom_dict[chrom][mirna][0] - 56 - flanks
                    end = chrom_dict[chrom][mirna][0] + flanks

                elif mirna_strand == '-':
                    start = chrom_dict[chrom][mirna][0] - flanks
                    end = chrom_dict[chrom][mirna][0] + 56 + flanks

                for read in alignments.fetch(chrom, start, end):
                    tally += 1
                    if tally % 1e5 == 0:
                        print('Lines read:',f'{tally:,}')
                    if read.has_tag("ts:i"):
                        if read.has_tag('UB') & read.has_tag('CB'):
                            count_table[mirna].add(read.get_tag("UB"))
                            MI_READS.write(read)
                            if read.cigarstring != '56M':
                                example_read = read

                            if read.is_reverse:
                                distance = read.reference_end - DROSHA_SITE
                            else:
                                distance = read.reference_start - DROSHA_SITE

                            # write detailed_results dictionary
                            if mirna not in detailed_results.keys():
                                detailed_results[mirna] = collections.defaultdict(list)

                            read_name = read.query_name
                            read_details = read.to_dict()
                            read_details.pop('name')
                            read_details.pop('tags')

                            detailed_results[mirna][read_name] = read_details

                            for tag in read.tags:
                                detailed_results[mirna][read_name][tag[0]] = tag[1]

                            detailed_results[mirna][read_name]['CB_UB'] = read.get_tag('CB') + '_' + read.get_tag('UB')
                            detailed_results[mirna][read_name]['Dist_to_DROSHA'] = distance

                        else:
                            NO_UB += 1
    MI_READS.close()
    alignments.close()


    ## Convert detailed results to a pandas dataframe
    print("Starting pandas dataframe conversion: ", time.asctime())

    results = pd.DataFrame()
    for m in detailed_results.keys():
        tmpdf = pd.DataFrame.from_dict(detailed_results[m], orient='index')
        tmpdf['miRNA'] = m
        results = pd.concat([tmpdf,results])
    print("Finish pandas dataframe conversion: ", time.asctime())

    print("Starting .bai index file writing: ", time.asctime())
    ## Write the .bai index file
    if os.path.exists(OUTBAM_FILENAME):

        print('Sorting output...')
        pysam.sort("-o", OUTBAM_FILENAME, "tmp_miRNA_matching.bam")
        os.remove("tmp_miRNA_matching.bam")
        if not os.path.exists('tmp_miRNA_matching.bam'):
            print('Cleaning up temp files...')
        print('Generating BAM index...')
        print('File size = ',np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2,2),'MB')
        pysam.index(OUTBAM_FILENAME)
        MI_READS.close()
    print("Finish .bai index file writing: ", time.asctime())

    print('Counting UMIs...', time.asctime())
    for mir in count_table.keys():
        count_table[mir] = len(count_table[mir])

    print('# candidate reads with no UB tag:',NO_UB)
    count_table = pd.DataFrame(count_table, index=['count']).T

    print("Finish count_miRNAs: ", time.asctime())
    return count_table, example_read, results

####### STOLEN AND BORKED, STILL TODO
"""

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
"""

def miRNA_to_featureMatrix(count_miRNAs_result, raw_feature_bc_matrix):
    """
    Add miRNAs to raw feature bc matrix.

    Args:
        count_miRNAs_result (pandas df): Dataframe output from count_miRNAs() containing
            all the information from pysam read alignments.
        raw_feature_bc_matrix (anndata): Raw matrix produced by 10x Genomics' CellRanger.

    Returns:
        AnnData: Raw matrix with miRNAs tacked onto the variable columns (gene names)

    """
    # # https://stackoverflow.com/questions/22412033/python-pandas-pivot-table-count-frequency-in-one-column

    barcode_miRNA_df = count_miRNAs_result[['CB','miRNA']].pivot_table(index='CB', columns='miRNA', aggfunc=len, fill_value=0)
    print('miRNAs by UMI count')
    print(barcode_miRNA_df.sum(0).sort_values(ascending=False))

    # cast barcode_miRNA_df as anndata and match up the .var attributes
    barcode_miRNA_adata = sc.AnnData(barcode_miRNA_df)
    barcode_miRNA_adata.var['gene_ids'] = barcode_miRNA_adata.var_names
    barcode_miRNA_adata.var['feature_types'] = raw_feature_bc_matrix.var['feature_types'][0]
    barcode_miRNA_adata.var['genome'] = raw_feature_bc_matrix.var['genome'][0]
    # rename all mirbase genes for clarity
    barcode_miRNA_adata.var_names = list(micro for micro in barcode_miRNA_adata.var_names + '_DroshaProd')

    # transpose and join outer
    print('raw_feature_bc_matrix dimensions: ', raw_feature_bc_matrix)
    print('barcode_miRNA_adata dimensions: ', barcode_miRNA_adata)
    tmp_combine = raw_feature_bc_matrix.T.concatenate(barcode_miRNA_adata.T, join='outer', index_unique=None)
    print('transpose merged dimensions: ', tmp_combine)
    raw_with_miRNAs = tmp_combine.T
    print('raw_with_miRNAs dimensions: ', raw_with_miRNAs)
    print('raw_with_miRNAs var_names: ', raw_with_miRNAs.var_names)
    del tmp_combine
    return raw_with_miRNAs

def main(args):
    outdir = args.outdir
    genome = args.genome

    if os.path.isdir(outdir):
        #overwrite = input('\nOutput directory already exists. Overwrite? Y/N ')
        #if overwrite.lower() == 'n':
        #    exit(0)
        #elif overwrite.lower() == 'y':
        #    shutil.rmtree(outdir)
        # Commented out above because if used as script (as it is now), there is no user input
        print('Output directory already exists')
    os.mkdir(outdir)

    print('Reading in gff3 file ...')
    miRNA_anno_df = attributes_to_columns_from_mirbase_GFF3(args.REFERENCE_FILE)
    PARENT_DICT = dict_of_parent_names(miRNA_anno_df)
    chrom_dict = make_dictionary_of_miRNA_Drosha_coords(miRNA_anno_df, PARENT_DICT)

    print('Detecting if BAM user \'chr\' prefix...')
    chrom_dict = fix_chr_chromnames(chrom_dict, args.BAMFILE)
    count_table, example_read, results = count_miRNAs(args.BAMFILE, chrom_dict)
    print('Done with part 1!')

    raw_feature_bc_matrix = sc.read_10x_h5(args.raw)
    raw_feature_bc_matrix.var_names_make_unique()
    raw_with_miRNAs = miRNA_to_featureMatrix(results, raw_feature_bc_matrix)

    outfile=os.path.join(outdir,'raw_feature_matrix_with_miRNAs.h5ad')
    raw_with_miRNAs.write(outfile)

    # testing purposes
    foo = sc.read(outfile)
    print('miRNA table dimensions: ', foo)

    if args.results_table == True:
        print('Writing miRNAs_result_table')
        results.to_csv(os.path.join(outdir, 'miRNAs_result_table.csv'), index=True)

    print('Done with part 2!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('REFERENCE_FILE', help='miRBase gff3 file')
    parser.add_argument('--outdir', dest='outdir', help='Output folder', default="raw_feature_matrix_with_miRNAs")
    parser.add_argument('--genome', dest='genome', help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
    #need to implement .mtx
    #parser.add_argument('--mtx', dest='mtx', help='Write output in 10X mtx format', default=False)
    parser.add_argument('--raw', dest='raw', required=True, help='10x Genomics raw feature bc matrix to concatenate miRNAs to')
    parser.add_argument('--results_table', dest='results_table', help='Write out results table of reading miRNAs as csv', default=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    main(args)
