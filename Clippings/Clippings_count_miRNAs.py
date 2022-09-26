#!/usr/bin/env python
# $ -cwd
# $ -v PATH,LD_LIBRARY_PATH
"""
Clippings_count_miRNAs.py

Count mapped reads harboring TSO sub-subsequence. Only works with CellRanger 4.0+ outputs.


Inputs:
    BAMFILE = BAM file produced by Cellranger v4+ with TSO reads tagged with ts:i
    REFERENCE_FILE = gff3 file with miRNA coordinates downloaded from miRBase
    genome
    outdir
    raw = Raw matrix to tack miRNAs onto the variable columns
    # TODO: Add ability to use miRNA GTF/coordinate files produced elsewhere
    # TODO: Change behavior of miRNA GFF3 references: Instead, specify genome version and pull from a pre-packaged file
    Usage: Clippings_count_miRNA.py BAMFILE REFERENCE_FILE --genome --outdir --raw

Author: jpreall@cshl.edu
Date: 2021-07
"""

import sys
import argparse
import os
import time
import json
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
import anndata
import scanpy as sc

def read_mirbase_gff3(file):
    """
    Reads a GFF3-formatted file from miRBase into a pandas dataframe.
    Note, a specific column format is expected:
        1. seqname  --> chromosome name
        2. source --> unused data source column
        3. feature --> indicates mature miRNA or precursor
        4. start --> start genomic coordinate
        5. end --> end genomic coordinate
        6. score --> unused
        7. strand --> +/- genomic orientation
        8. frame --> unused (non-coding!)
        9. attributes --> ';' delimited additional attributes. Extensible.

        The 'attributes' column will be expanded to create N additional columns
        to capture all non-redundant attributes present in the GFF3 file.

    Args:
        file (string path): Path to .gff3 formatted file with miRNA annotations downloaded from miRBase.

    Returns:
        pandas DataFrame: Dataframe with columns containing the attributes from GFF3 file.

    """
    columns = [
        'seqname', 'source', 'feature', 'start', 
        'end', 'score', 'strand', 'frame', 'attributes']

    df = pd.read_table(file, comment='#', header=None)
    df.columns = columns

    split_attributes = df.attributes.apply(
        lambda attr_glob: dict(
            [attr.split(sep="=", maxsplit=1) for attr in attr_glob.split(";")]
        )
    )

    miRNA_anno_df = pd.concat(
        [
        df,
        pd.concat(
            pd.DataFrame(item, index=[idx]) 
            for idx, item in split_attributes.iteritems()
        ).fillna('None')
        ], 
        axis=1).drop(columns='attributes')

    return miRNA_anno_df

def make_Drosha_coord_dict(miRNA_anno_df):
    """
    Clippings will detect any TSO read within a few bases of the predicted Drosha cleavage site
    This is defined as the 3'-end of the 3p arm of the mature miRNA

    Args:
        miRNA_anno_df (pandas DataFrame): produced by read_mirbase_gff3()
        #parent_dict (dict): Dictionary of ID to Name of only primary transcripts.

    Returns:
        dict: sorted by chromosome of Drosha cleavage sites of all mature miRNAs.
        Default dict has miRNA name and drosha coordinate with '+' or '-' strand.
    """
    
    try:
        threep = miRNA_anno_df[miRNA_anno_df['Name'].str.match('.*3p$')]
    except KeyError:
        threep = miRNA_anno_df[miRNA_anno_df['gene_name'].str.match('.*3p$')]
    
    # Extract only 3p-arms of miRNAs
    threep = miRNA_anno_df[miRNA_anno_df['Name'].str.match('.*3p$', re.IGNORECASE)].copy()
    
    # De-duplicate any miRNAs with the same mature name, if they come from different parents
    for idx, count in threep.groupby('Name').cumcount().iteritems():
        if count != 0:
            threep.loc[idx,'Name'] = threep.loc[idx,'Name'].replace('-3p',f'.{count}-3p')

    def Drosha_site(gff_3p_row):
        COORD = gff_3p_row['start'] if \
        gff_3p_row['strand'] == '-' \
        else gff_3p_row['end']
        
        return (COORD, gff_3p_row['strand'])

    coord_dict = defaultdict(dict)
    for i, row in threep.iterrows():
        coord_dict[row.seqname][row.Name] = Drosha_site(row)
    
    return coord_dict
    
def fix_chromnames(coord_dict, BAM):
    """
    Checks if BAM file is using 'chr' prefix on chromosome names and fixes the dictionary if necessary.

    Args:
        coord_dict (dict): Dictionary of chromosome as key and value as default dict.
            Default dict has miRNA name and drosha coordinate with '+' or '-' strand.
        BAM (bamfile): Bamfile produced by 10x Genomics' CellRanger

    Returns:
        dict: Updated version of coord_dict with correct 'chr' prefix.

    """
    chromnames_in_gff3 = list(coord_dict.keys())
    chromnames_in_bam = pysam.AlignmentFile(BAM, "rb").header.references
    pattern = re.compile('chr', re.IGNORECASE)
    
    def get_prefix(list_of_chromosome_names):
        dominant_prefix = pd.Series([
            re.match(pattern,name).group() \
            for name in list_of_chromosome_names \
            if re.match(pattern, name)
        ]
        ).value_counts().idxmax()
        
        return dominant_prefix
    
    # Check to see if either the BAM file or the GFF using a 'chr'-style prefix
    BAM_has_prefix = pd.Series(chromnames_in_bam).str.match(pattern).sum() / len(chromnames_in_bam) > .1
    BAM_prefix = get_prefix(chromnames_in_bam) if BAM_has_prefix else ''
    
    GFF_has_prefix = pd.Series(chromnames_in_gff3).str.match(pattern).sum() / len(chromnames_in_bam) > .1
    GFF_prefix = get_prefix(chromnames_in_gff3) if GFF_has_prefix else ''

    ## rewrite GFF prefix to match BAM
    if BAM_has_prefix and not GFF_has_prefix:
        for c in chromnames_in_gff3:
            fixed_c = BAM_prefix + str(c)
            if fixed_c in chromnames_in_bam:
                print(f'Changing {c} to {fixed_c}')
                coord_dict[fixed_c] = coord_dict.pop(c)
                
    elif GFF_has_prefix and not BAM_has_prefix:
        for c in chromnames_in_gff3:
            fixed_c = c.replace(GFF_prefix,BAM_prefix)
            if fixed_c in chromnames_in_bam:
                print(f'Changing {c} in GFF reference to {fixed_c}')
                coord_dict[fixed_c] = coord_dict.pop(c)
    else:
        print('BAM and GFF prefixes match. No fixing needed.')
            
    return coord_dict

def get_sample_name_from_bam_header(BAMFILE):
    alignments = pysam.AlignmentFile(BAMFILE, "rb")    
    try: 
        sample_name = alignments.header.get('RG')[0].get('SM')
    except:
        sample_name = 'UnknownSample'
    return sample_name

def get_bam_readlength(BAMFILE):
    alignments = pysam.AlignmentFile(BAMFILE, 'rb')
    rlens = []
    for n, read in enumerate(alignments.fetch()):
        rlens += [len(read.get_forward_sequence())]
        if n > 10000:
            break
    return pd.Series(rlens).value_counts().idxmax()

def get_bam_features(BAMFILE):
    # Get Sample name from Cellranger-formatted BAM:
    alignments = pysam.AlignmentFile(BAMFILE, "rb")
    try: 
        sample_name = alignments.header.get('RG')[0].get('SM')
    except:
        sample_name = 'UnknownSample'
    
    # Get read length of run by skimming first 10,000 lines:
    rlens = []
    for n, read in enumerate(alignments.fetch()):
        rlens += [len(read.get_forward_sequence())]
        if n > 10000:
            break
    read_length = pd.Series(rlens).value_counts().idxmax()
    
    return sample_name, read_length
    
def count_miRNAs(BAM, chrom_dict, sampleName, flanks=0):
    """
    Reads bamfile and stores read information into a pandas DataFrame and counts number
    of each miRNA. Also outputs miRNA labelled reads as a bam file with index.

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
            miRNA_dist_dict (dict): Dictionary with distance to Drosha counts for each miRNA.

    """
    print("Starting count_miRNAs: ", time.asctime())
    example_read = None

    allmi = {}
    for c in chrom_dict.keys():
        for m in chrom_dict[c].keys():
            allmi[m] = c

    print(BAM)
    alignments = pysam.AlignmentFile(BAM, "rb")
    MI_READS = pysam.AlignmentFile('tmp_miRNA_matching.bam', "wb", template=alignments)

    count_table = defaultdict(set)
    detailed_results = defaultdict(list)
    miRNA_dist_dict = defaultdict(dict)

    tally = 0
    NO_UB = 0
    bam_read_length = 56

    Drosha_tolerance = 4
    dist_DROSHA_dict = dict.fromkeys(range(-Drosha_tolerance, Drosha_tolerance+1), 0)

    for mirna, chrom in allmi.items():
        chroms = [chrom]
        mirnas = [mirna]
       
        miRNA_dist_dict[mirna] = dist_DROSHA_dict.copy()

        mirna_strand = chrom_dict[chrom][mirna][1]
        DROSHA_SITE = chrom_dict[chrom][mirna][0]

        if mirna_strand == '-':
            start = DROSHA_SITE - bam_read_length - flanks
            end = DROSHA_SITE + flanks

        elif mirna_strand == '+':
            start = DROSHA_SITE - flanks
            end = DROSHA_SITE + bam_read_length + flanks

        for read in alignments.fetch(chrom, start, end):
            tally += 1
            if tally % 1e5 == 0:
                print('Lines read:', f'{tally:,}')
                # print('count_table: ', count_table)

            if read.has_tag("ts:i"):
                if read.has_tag('UB') & read.has_tag('CB'):
                    if read.is_reverse:  # basically '-' strand
                        #distance = read.reference_end - DROSHA_SITE
                        distance = DROSHA_SITE - read.reference_end
                    else:  # basically '+' strand
                        distance = read.reference_start - DROSHA_SITE
                    # distance to DROSHA should be no less than 'X'
                    if abs(distance) < 5:
                        initial_len = len(count_table[mirna])
                        count_table[mirna].add(read.get_tag("UB"))
                        # counting the UB's distance to DROSHA if condition met
                        if len(count_table[mirna]) != initial_len:
                            miRNA_dist_dict[mirna][distance] += 1

                        MI_READS.write(read)

                        if read.cigarstring != '56M':
                            example_read = read

                        # write detailed_results dictionary
                        if mirna not in detailed_results.keys():
                            detailed_results[mirna] = defaultdict(list)

                        read_name = read.query_name
                        read_details = read.to_dict()
                        read_details.pop('name')
                        read_details.pop('tags')

                        detailed_results[mirna][read_name] = read_details

                        for tag in read.tags:
                            detailed_results[mirna][read_name][tag[0]] = tag[1]

                        detailed_results[mirna][read_name]['CB_UB'] = read.get_tag(
                            'CB') + '_' + read.get_tag('UB')
                        detailed_results[mirna][read_name]['Dist_to_DROSHA'] = distance

                else:
                    NO_UB += 1

    MI_READS.close()
    alignments.close()

    # Convert detailed results to a pandas dataframe
    print("Starting pandas dataframe conversion: ", time.asctime())

    results = pd.DataFrame()
    print(type(detailed_results))
    for key, item in detailed_results.items():
        tmpdf = pd.DataFrame.from_dict(item, orient='index')
        tmpdf['miRNA'] = key
        results = pd.concat([tmpdf,results])

    #for m in detailed_results.keys():
        #tmpdf = pd.DataFrame.from_dict(detailed_results[m], orient='index')
        #tmpdf['miRNA'] = m
        #results = pd.concat([tmpdf, results])
    print("Finish pandas dataframe conversion: ", time.asctime())

    print("Starting .bai index file writing: ", time.asctime())
    # Write the .bai index file
    print("Sorting output...")
    OUTBAM_FILENAME = sampleName + '_sorted_miRNA_reads.bam'
    pysam.sort("-o", OUTBAM_FILENAME, "tmp_miRNA_matching.bam")
    os.remove("tmp_miRNA_matching.bam")
    if os.path.exists('tmp_miRNA_matching.bam'):
        os.remove("tmp_miRNA_matching.bam")
        print('Cleaning up temp files...')
    print('Generating BAM index...', time.asctime())
    print('File size = ', np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2, 2), 'MB')
    pysam.index(OUTBAM_FILENAME)
    MI_READS.close()
    print("Finish .bai index file writing: ", time.asctime())

    print('Counting UMIs...', time.asctime())
    for mir in count_table.keys():
        count_table[mir] = len(count_table[mir])

    print('# candidate reads with no UB tag:', NO_UB)
    count_table = pd.DataFrame(count_table, index=['count']).T

    print("Finish count_miRNAs: ", time.asctime())
    return count_table, example_read, results, miRNA_dist_dict



def count_miRNAs_OLD(BAM, chrom_dict, sampleName, flanks=0):
    """
    Reads bamfile and stores read information into a pandas DataFrame and counts number
    of each miRNA. Also outputs miRNA labelled reads as a bam file with index.

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
            miRNA_dist_dict (dict): Dictionary with distance to Drosha counts for each miRNA.

    """
    print("Starting count_miRNAs: ", time.asctime())
    example_read = None

    allmi = {}
    for c in chrom_dict.keys():
        for m in chrom_dict[c].keys():
            allmi[m] = c

    print(BAM)
    alignments = pysam.AlignmentFile(BAM, "rb")
    MI_READS = pysam.AlignmentFile('tmp_miRNA_matching.bam', "wb", template=alignments)

    count_table = defaultdict(set)
    detailed_results = defaultdict(list)
    tally = 0
    NO_UB = 0

    miRNA_dist_dict = defaultdict(dict)

    for i in range(len(allmi.keys())):
        m = list(allmi.keys())[i]
        chroms = [allmi[m]]
        mirnas = [m]

        dist_DROSHA_dict = dict.fromkeys(range(-4, 5), 0)
        miRNA_dist_dict[m] = dist_DROSHA_dict
        #top_miRNAs_sc = ['hsa-mir-16-2', 'hsa-mir-10394', 'hsa-mir-4709', 'hsa-mir-6805', 'hsa-mir-10393', 'hsa-mir-6797',
        #                 'hsa-mir-103a-2', 'hsa-mir-7111', 'hsa-mir-101-1', 'hsa-mir-21']

        #if m in top_miRNAs_sc:
        #    print('top_miRNAs_sc: ', m)

        for chrom in chroms:
            for mirna in mirnas:

                mirna_strand = chrom_dict[chrom][mirna][1]
                DROSHA_SITE = chrom_dict[chrom][mirna][0]

                if mirna_strand == '-':
                    start = chrom_dict[chrom][mirna][0] - 56 - flanks
                    end = chrom_dict[chrom][mirna][0] + flanks

                elif mirna_strand == '+':
                    start = chrom_dict[chrom][mirna][0] - flanks
                    end = chrom_dict[chrom][mirna][0] + 56 + flanks

                for read in alignments.fetch(chrom, start, end):
                    tally += 1
                    if tally % 1e5 == 0:
                        print('Lines read:', f'{tally:,}')
                        # print('count_table: ', count_table)

                    if read.has_tag("ts:i"):
                        if read.has_tag('UB') & read.has_tag('CB'):
                            if read.is_reverse:  # basically '-' strand
                                #distance = read.reference_end - DROSHA_SITE
                                distance = DROSHA_SITE - read.reference_end
                            else:  # basically '+' strand
                                distance = read.reference_start - DROSHA_SITE
                            # distance to DROSHA should be no less than 'X'
                            if abs(distance) < 5:
                                initial_len = len(count_table[mirna])
                                count_table[mirna].add(read.get_tag("UB"))
                                # counting the UB's distance to DROSHA if condition met
                                if len(count_table[mirna]) != initial_len:
                                    miRNA_dist_dict[m][distance] += 1

                                MI_READS.write(read)

                                if read.cigarstring != '56M':
                                    example_read = read

                                # write detailed_results dictionary
                                if mirna not in detailed_results.keys():
                                    detailed_results[mirna] = defaultdict(list)

                                read_name = read.query_name
                                read_details = read.to_dict()
                                read_details.pop('name')
                                read_details.pop('tags')

                                detailed_results[mirna][read_name] = read_details

                                for tag in read.tags:
                                    detailed_results[mirna][read_name][tag[0]] = tag[1]

                                detailed_results[mirna][read_name]['CB_UB'] = read.get_tag(
                                    'CB') + '_' + read.get_tag('UB')
                                detailed_results[mirna][read_name]['Dist_to_DROSHA'] = distance

                        else:
                            NO_UB += 1

    MI_READS.close()
    alignments.close()

    # Convert detailed results to a pandas dataframe
    print("Starting pandas dataframe conversion: ", time.asctime())

    results = pd.DataFrame()
    print(type(detailed_results))
    for key, item in detailed_results.items():
        tmpdf = pd.DataFrame.from_dict(item, orient='index')
        tmpdf['miRNA'] = key
        results = pd.concat([tmpdf,results])

    #for m in detailed_results.keys():
        #tmpdf = pd.DataFrame.from_dict(detailed_results[m], orient='index')
        #tmpdf['miRNA'] = m
        #results = pd.concat([tmpdf, results])
    print("Finish pandas dataframe conversion: ", time.asctime())

    print("Starting .bai index file writing: ", time.asctime())
    # Write the .bai index file
    print("Sorting output...")
    OUTBAM_FILENAME = sampleName + '_sorted_miRNA_reads.bam'
    pysam.sort("-o", OUTBAM_FILENAME, "tmp_miRNA_matching.bam")
    os.remove("tmp_miRNA_matching.bam")
    if os.path.exists('tmp_miRNA_matching.bam'):
        os.remove("tmp_miRNA_matching.bam")
        print('Cleaning up temp files...')
    print('Generating BAM index...', time.asctime())
    print('File size = ', np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2, 2), 'MB')
    pysam.index(OUTBAM_FILENAME)
    MI_READS.close()
    print("Finish .bai index file writing: ", time.asctime())

    print('Counting UMIs...', time.asctime())
    for mir in count_table.keys():
        count_table[mir] = len(count_table[mir])

    print('# candidate reads with no UB tag:', NO_UB)
    count_table = pd.DataFrame(count_table, index=['count']).T

    print("Finish count_miRNAs: ", time.asctime())
    return count_table, example_read, results, miRNA_dist_dict


def dist_toDrosha_table(path_to_json, min_counts=10, scaleby='total'):
    """
    Creates distance to drosha count table for miRNAs scaled by max value or total counts for given a json file,
    which is then used downstream to label bogus miRNAs.

    Args:
        path_to_json (str): Path to json file with dictionary of distance to Drosha counts for each miRNA.
        min_counts (int): Min counts to consider for labeling bogus miRNAs. Must be at least 5 to calculate pileups.
        scaleby (str): Method to scale counts to 0-1.

    Returns:
        droshaDist_sliced_norm (Pandas Dataframe): Dataframe with scaled counts of distance to drosha for each miRNA.
        To be used for labeling bogus miRNAs.

    """
    if min_counts < 5:
        raise ValueError('min_counts must be at least 5 in order to label bogus miRNAs')
    def data_loader(path_to_json):
        f = open(path_to_json)
        # returns JSON object as
        # a dictionary
        data = json.load(f)
        f.close()
        droshaDist = pd.DataFrame(data).T
        print('top 10 miRNA counts', droshaDist.sum(1).sort_values(ascending=False).head(10))
        return droshaDist

    droshaDist = data_loader(path_to_json)
    droshaDist_sliced = droshaDist[droshaDist.sum(1) >= min_counts].copy()
    droshaDist_sliced.loc[:, 'rowSum'] = droshaDist_sliced.sum(1)
    droshaDist_sliced = droshaDist_sliced.sort_values(by='rowSum', ascending=False)
    del droshaDist_sliced['rowSum']

    if scaleby == 'max':
        droshaDist_sliced_norm = droshaDist_sliced.div(droshaDist_sliced.max(axis=1), axis=0)
    elif scaleby == 'total':
        droshaDist_sliced_norm = droshaDist_sliced.div(droshaDist_sliced.sum(axis=1), axis=0)
    elif scaleby == 'none':
        droshaDist_sliced_norm = droshaDist
    else:
        raise ValueError('scaleby must be one of the following: max, total, none')
    # fill in any na values
    droshaDist_sliced_norm = droshaDist_sliced_norm.fillna(0)

    #sns.clustermap(droshaDist_sliced_norm, col_cluster=False, row_cluster=False, yticklabels=True)
    return droshaDist_sliced_norm


def largest_index(arr):
    """
    Find the index with the largest value in an array

    Args:
        arr (array_like): Array containing number with largest value desired.

    Returns:
       index_of_max (int): the index of the largest number in arr.

    """
    if len(arr) < 9:
        raise ValueError('length of array must be greater than 9 in order to find DROSHA pileup pattern')
    midpoint = len(arr) // 2
    # only look at center near DROSHA site
    index_of_max = 0
    for i in range(1, len(arr), 1):
        if arr[i] > arr[index_of_max]:
            index_of_max = i
    return index_of_max


def label_bogus(arr):
    """
       Label bogus miRNAs based on Drosha site read pileups.

       Args:
           arr (array_like): Array containing relative drosha site pileup information.

       Returns:
           boolean: True/False if the miRNA should be considered bogus.

    """
    peak = largest_index(arr)
    if arr[peak] < 0.5:
        return True  # if largest pileup < 50% then miRNA is bogus
    else:
        if sum(arr[peak:peak + 3]) > 0.75:
            return False  # if the right side (up to 2bp away) DROSHA site pileup is greater than 75% then not bogus
        elif sum(arr[peak - 2:peak]) > 0.75:
            return False  # if the left side (up to 2bp away) DROSHA site pileup is greater than 75% then not bogus
        else:
            return True  # if no distinct pileup, then miRNA is considered bogus


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
    # remove duplicated UBs from the result table that has all the queries
    count_miRNAs_result = count_miRNAs_result.drop_duplicates(subset=['UB'])
    barcode_miRNA_df = count_miRNAs_result[['CB', 'miRNA']].pivot_table(
        index='CB', columns='miRNA', aggfunc=len, fill_value=0)
    print('miRNAs by UMI count')
    print(barcode_miRNA_df.sum(0).sort_values(ascending=False).head(10))

    # cast barcode_miRNA_df as anndata and match up the .var attributes
    barcode_miRNA_adata = sc.AnnData(barcode_miRNA_df)
    barcode_miRNA_adata.var['gene_ids'] = barcode_miRNA_adata.var_names
    barcode_miRNA_adata.var['feature_types'] = raw_feature_bc_matrix.var['feature_types'][0]
    barcode_miRNA_adata.var['genome'] = raw_feature_bc_matrix.var['genome'][0]
    # rename all mirbase genes for clarity
    barcode_miRNA_adata.var_names = list(
        micro for micro in barcode_miRNA_adata.var_names + '_DroshaProd')

    # transpose and join outer
    print('raw_feature_bc_matrix dimensions: ', raw_feature_bc_matrix)
    print('barcode_miRNA_adata dimensions: ', barcode_miRNA_adata)
    tmp_combine = raw_feature_bc_matrix.T.concatenate(
        barcode_miRNA_adata.T, join='outer', index_unique=None)
    print('transpose merged dimensions: ', tmp_combine)
    raw_with_miRNAs = tmp_combine.T
    print('raw_with_miRNAs dimensions: ', raw_with_miRNAs)
    print('raw_with_miRNAs var_names: ', raw_with_miRNAs.var_names)
    del tmp_combine
    return raw_with_miRNAs


def main(cmdl):
    args = _parse_cmdl(cmdl)
    outdir = args.outdir
    genome = args.genome
    if os.path.isdir(outdir):
        # overwrite = input('\nOutput directory already exists. Overwrite? Y/N ')
        # if overwrite.lower() == 'n':
        #    exit(0)
        # elif overwrite.lower() == 'y':
        #    shutil.rmtree(outdir)
        # Commented out above because if used as script (as it is now), there is no user input
        print('Output directory already exists')
    os.mkdir(outdir)
    
    sampleName = get_sample_name_from_bam_header(args.BAMFILE)
    print(f'Detected sample name as: {sampleName}')

    print('Reading in gff3 file ...')
    miRNA_anno_df = read_mirbase_gff3(args.REFERENCE_FILE)
    coord_dict = make_Drosha_coord_dict(miRNA_anno_df)


    print('Detecting if BAM uses \'chr\' prefix...')
    coord_dict = fix_chromnames(coord_dict, args.BAMFILE)
    count_table, example_read, results, miRNA_dist_dict = count_miRNAs(args.BAMFILE, coord_dict, sampleName)
    with open(os.path.join(outdir, "distance_to_DROSHA_dict.json"), "w") as outfile:
        json.dump(miRNA_dist_dict, outfile)
    print('Done with part 1!')

    raw_feature_bc_matrix = sc.read_10x_h5(args.raw)
    raw_feature_bc_matrix.var_names_make_unique()
    raw_with_miRNAs = miRNA_to_featureMatrix(results, raw_feature_bc_matrix)
    raw_with_miRNAs.var['sample_name'] = sampleName

    # label bogus miRNAs in drosha dist table
    drosha_table = dist_toDrosha_table(os.path.join(outdir, "distance_to_DROSHA_dict.json"))
    drosha_table.index += '_DroshaProd'  # rename index so we can merge table with raw_with_miRNAs
    drosha_table['bogus_miRNA'] = drosha_table.apply(lambda row: label_bogus(row), axis=1)

    # merge drosha dist table with raw_with_miRNAs table
    raw_with_miRNAs.var = raw_with_miRNAs.var.join(drosha_table['bogus_miRNA'])

    outfile = os.path.join(outdir, 'raw_feature_matrix_with_miRNAs.h5ad')
    raw_with_miRNAs.write(outfile)

    # testing purposes
    #foo = sc.read(outfile)
    #print('miRNA table dimensions: ', foo)

    if args.results_table is True:
        print('Writing miRNAs_result_table')
        results.to_csv(os.path.join(outdir, 'miRNAs_result_table.csv'), index=True)

    print('Done with part 2!')

def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Clippings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('REFERENCE_FILE', help='miRBase gff3 file')
    parser.add_argument('--outdir', dest='outdir', help='Output folder',
                        default="raw_feature_matrix_with_miRNAs")
    parser.add_argument('--genome', dest='genome',
                        help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
    # need to implement .mtx
    # parser.add_argument('--mtx', dest='mtx', help='Write output in 10X mtx format', default=False)
    parser.add_argument('--raw', dest='raw', required=True,
                        help='10x Genomics raw feature bc matrix to concatenate miRNAs to')
    parser.add_argument('--results_table', dest='results_table',
                        help='Write out results table of reading miRNAs as csv', default=True)
    #parser = logmuse.add_logging_options(parser)
    return parser.parse_args(cmdl)

if __name__ == '__main__':
    # Override sys.argv
#    sys.argv = ['Clippings_count_miRNAs.py', '/mnt/grid/scc/data/Preall/Preall_CR01/count/Preall_CR01_H_neg/outs/possorted_genome_bam.bam',
#                '/grid/preall/home/bhe/microRNA_project/hsa.gff3', '--outdir', 'testing_mar22_CR01_H_neg',
#                '--genome', '/mnt/grid/scc/data/CellRanger/references/refdata-gex-GRCh38-2020-A/',
#                '--raw', '/mnt/grid/scc/data/Preall/Preall_CR01/count/Preall_CR01_S_plus/outs/raw_feature_bc_matrix.h5']
    print(sys.argv[1:])
    main(sys.argv[1:])
