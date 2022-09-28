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
from typing import Union

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
    alignments.close()
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
    bam_read_length = pd.Series(rlens).value_counts().idxmax()
    
    return sample_name, bam_read_length
    
def count_miRNAs(BAM, coord_dict, sampleName, bam_read_length, write_miRNA_bam=False):
    """
    Reads bamfile and stores read information into a pandas DataFrame and counts number
    of each miRNA. Also outputs miRNA labelled reads as a bam file with index.

    Args:
        BAM (bamfile): Bamfile produced by 10x Genomics' CellRanger
        
        coord_dict (dict): Dictionary of miRNA cropping loci sorted by chromosome. 
            Produced by make_Drosha_coord_dict()

    Returns:
        (tuple): tuple containing:

            count_table (pandas DataFrame): Dataframe of counts of each miRNA.
            example_read (pysam AlignmentSegment): Read as AlignmentSegment.
            results (pandas DataFrame): Dataframe of all reads and info from pysam.
            miRNA_dist_dict (dict): Dictionary with distance to Drosha counts for each miRNA.

    """
    print("Starting count_miRNAs: ", time.asctime())

    allmi = {}
    for c in coord_dict.keys():
        for m in coord_dict[c].keys():
            allmi[m] = c

    print(BAM)
    alignments = pysam.AlignmentFile(BAM, "rb")
    
    if write_miRNA_bam:
        MI_READS = pysam.AlignmentFile('tmp_miRNA_matching.bam', "wb", template=alignments)

    count_table = defaultdict(set)
    detailed_results = defaultdict(list)
    miRNA_dist_dict = defaultdict(dict)
    DEBUG_read_details = defaultdict(dict)
    
    NO_UB = 0

    Drosha_tolerance = 4
    dist_DROSHA_dict = dict.fromkeys(range(-Drosha_tolerance, Drosha_tolerance+1), 0)

    for mirna, chrom in allmi.items():
       
        miRNA_dist_dict[mirna] = dist_DROSHA_dict.copy()

        mirna_strand = coord_dict[chrom][mirna][1]
        DROSHA_SITE = coord_dict[chrom][mirna][0]

        if mirna_strand == '-':
            start = DROSHA_SITE - bam_read_length
            end = DROSHA_SITE 

        elif mirna_strand == '+':
            start = DROSHA_SITE
            end = DROSHA_SITE + bam_read_length
            
        
        # Inspect all reads in the vicinity of the expected Drosha cropping site
        for tally, read in enumerate(alignments.fetch(chrom, start, end)):
            tally += 1
            if tally % 1e5 == 0:
                print('Lines read:', f'{tally:,}')

            if read.has_tag("ts:i"):
                if read.has_tag('UB') & read.has_tag('CB'):
                    
                    # Compute distance from TSO incorporation site to predicted Drosha cropping site
                    if read.is_reverse:  # '-' strand
                        #distance = read.reference_end - DROSHA_SITE
                        distance = DROSHA_SITE - read.reference_end
                    else:  # '+' strand
                        distance = read.reference_start - DROSHA_SITE
                        
                    # Distance to cropping site should be no less than 'Drosha_tolerance'
                    if abs(distance) <= Drosha_tolerance:
                        UB = read.get_tag("UB")
                        
                        # keep a running tally of distances to cropping site per unique miRNA count
                        if UB not in count_table[mirna]:
                            miRNA_dist_dict[mirna][distance] += 1
                        
                        # Add this UMI to the the set
                        count_table[mirna].add(UB)
                        
                        # Write matching reads to BAM output
                        if write_miRNA_bam:
                            MI_READS.write(read)

                        # write detailed_results dictionary
                        if mirna not in detailed_results.keys():
                            detailed_results[mirna] = defaultdict(list)

                        read_name = read.query_name
                        DEBUG_read_details[read_name] = read.to_dict()
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

    if write_miRNA_bam:
        MI_READS.close()
        
    alignments.close()

    # Convert detailed results to a pandas dataframe
    print("Starting pandas dataframe conversion: ", time.asctime())

    results = pd.DataFrame()
    for key, item in detailed_results.items():
        tmpdf = pd.DataFrame.from_dict(item, orient='index')
        tmpdf['miRNA'] = key
        results = pd.concat([tmpdf,results])

    print("Finish pandas dataframe conversion: ", time.asctime())
    
    if write_miRNA_bam:
        print("Finalizing BAM output: ", time.asctime())
        # Write the .bai index file
        
        print("Sorting output...")
        
        #Note: This file should already be sorted if coming from Cellranger.
        OUTBAM_FILENAME = sampleName + '_sorted_miRNA_reads.bam'
        pysam.sort("-o", OUTBAM_FILENAME, "tmp_miRNA_matching.bam")

        if os.path.exists('tmp_miRNA_matching.bam'):
            os.remove("tmp_miRNA_matching.bam")
            print('Cleaning up temp files...')
            
        FILESIZE = np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2, 2)
        print(f'File size = {FILESIZE} MB')

        print('Generating BAM index...', time.asctime())
        pysam.index(OUTBAM_FILENAME)
        print("Finished sorting and indexing BAM output ", time.asctime())

    print('Counting total UMIs per miRNA...', time.asctime())
    #counts_per_miRNA = pd.Series({key:len(count_table[key]) for key in count_table.keys()})
    counts_per_miRNA = pd.DataFrame(
        {key:len(UMIlist) for key, UMIlist in count_table.items()},
        index=['counts']).T

    print('# candidate reads with no UB tag:', NO_UB)

    print("Finished count_miRNAs: ", time.asctime())
    #return count_table, results, miRNA_dist_dict, detailed_results, read, counts_per_miRNA, DEBUG_read_details
    return results, miRNA_dist_dict, counts_per_miRNA

def count_miRNAs_simple(BAM: Union[str, bytes, os.PathLike], 
                        coord_dict: dict, 
                        sample_name: str, 
                        bam_read_length: int,
                        Drosha_tolerance: int=4,
                        write_miRNA_bam: bool=False):
    """
    Reads bamfile and stores read information into a pandas DataFrame and counts number
    of each miRNA. Also outputs miRNA labelled reads as a bam file with index.

    Args:
        BAM (bamfile): BAM-formatted alignment file produced by 10x Genomics' 
             CellRanger pipeline to be scanned for miRNA processing products.
        
        coord_dict (dict): Dictionary of miRNA cropping loci sorted by chromosome. 
            Produced by make_Drosha_coord_dict()
            
        bam_read_length (int): Read length detected in BAM files. Used to extract
            reads within 1 read length of Drosha coordinates for counting.
            
        Drosha_tolerance: Maximum distance (inclusive) between predicted Drosha cleavage
            site and read TSO-incorporation site to include read in miRNA count
        
        write_miRNA_bam (bool): If True, writes ALL TSO-containing reads within a
             window of +/- one read length (detected from BAM file) to a BAM file

    Returns:
        (tuple): tuple containing:
            detailed_results_df (pd.DataFrame): All TSO-containing reads in a BAM-like
                dataframe, including nearest miRNA name and distance to Drosha site
            count_df (pd.DataFrame): Unique cell-by-gene count matrix for all miRNAs.
            counts_per_miRNA (pd.Series): Sorted total abundance of each miRNA detected
    """
    
    def stash_read(
        read: pysam.AlignedSegment, 
        detailed_results: dict,
        miRNA: str,
        distance: int):
        """
        Any read passing filter as a putative Drosha cropping product 
        is stashed in a nested dictionary
        
        Args:
            read: read object generated by pysam
            detailed_results (dict): dictionary in which results are stashed
            miRNA (str): Name of miRNA in vicinity of read
            distance (int): Distance from read 5'-end to predicted Drosha cropping site
        """
    
        read_name = read.query_name
        read_details = read.to_dict()
        read_details.pop('tags')
        read_details.pop('name')

        detailed_results[read_name] = read_details
        detailed_results[read_name]['miRNA'] = mirna
        detailed_results[read_name]['CB_UB'] = read.get_tag(
            'CB') + '_' + read.get_tag('UB')
        detailed_results[read_name]['Dist_to_DROSHA'] = distance

        for tag in read.tags:
            detailed_results[read_name][tag[0]] = tag[1]
            
    def finalize_miRNA_bam(
        bamfile: Union[str, bytes, os.PathLike] = 'tmp_miRNA_matching.bam',
        sample_name: str = sample_name):
        """
        Sorts, indexes, and cleans up temporary miRNA BAM file.
        """
        print("Finalizing BAM output: ", time.asctime())    
        print("Sorting output...")
        
        OUTBAM_FILENAME = sample_name + '_sorted_miRNA_reads.bam'
        pysam.sort("-o", OUTBAM_FILENAME, "tmp_miRNA_matching.bam")

        if os.path.exists('tmp_miRNA_matching.bam'):
            os.remove("tmp_miRNA_matching.bam")
            print('Cleaning up temp files...')
            
        FILESIZE = np.round(os.path.getsize(OUTBAM_FILENAME) / 1024**2, 2)
        print(f'File size = {FILESIZE} MB')

        print('Generating BAM index...', time.asctime())
        pysam.index(OUTBAM_FILENAME)
        print("Finished sorting and indexing BAM output ", time.asctime())
        
    print("Starting count_miRNAs: ", time.asctime())
    #count_table = defaultdict(set)
    #miRNA_dist_dict = defaultdict(dict)
    
    # Setup output variables
    detailed_results = defaultdict(list)
    dist_DROSHA_dict = dict.fromkeys(range(-Drosha_tolerance, Drosha_tolerance+1), 0)
    NO_UB = 0
    
    # Create dict of all miRNA names and their respective chromosomes to iterate over
    allmi = {}
    for c in coord_dict.keys():
        for m in coord_dict[c].keys():
            allmi[m] = c

    print(f'Reading BAM file: {BAM} ...') 
    print(f'Inferred read length: {bam_read_length} bases')
    alignments = pysam.AlignmentFile(BAM, "rb")
    
    # Initialize output BAM file, if writing one
    if write_miRNA_bam:
        MI_READS = pysam.AlignmentFile('tmp_miRNA_matching.bam', "wb", template=alignments)

    # Main loop over all miRNAs in GTF file to detect candidate fragments:
    for mirna, chrom in allmi.items():

        mirna_strand = coord_dict[chrom][mirna][1]
        DROSHA_SITE = coord_dict[chrom][mirna][0]

        if mirna_strand == '-':
            start = DROSHA_SITE - bam_read_length
            end = DROSHA_SITE 

        elif mirna_strand == '+':
            start = DROSHA_SITE
            end = DROSHA_SITE + bam_read_length


        # Inspect all reads in the vicinity of the expected Drosha cropping site
        for tally, read in enumerate(alignments.fetch(chrom, start, end)):
            tally += 1
            if tally % 1e5 == 0:
                print('Lines read:', f'{tally:,}')

            if read.has_tag("ts:i"):
                if read.has_tag('UB') & read.has_tag('CB'):

                    # Compute distance from TSO incorporation site to predicted Drosha cropping site
                    if read.is_reverse:  # '-' strand
                        distance = DROSHA_SITE - read.reference_end
                        
                    else:  # '+' strand
                        distance = read.reference_start - DROSHA_SITE
                    
                    #Stash this read in the results dictionary
                    stash_read(read, detailed_results, mirna, distance)
                    
                    # Write matching reads to BAM output
                    if write_miRNA_bam:
                        MI_READS.write(read)
                        
                else:
                    NO_UB += 1
    alignments.close()
    
    print('# candidate reads with no UB tag:', NO_UB)
    
    if write_miRNA_bam:
        MI_READS.close()
        finalize_miRNA_bam()
        
    # Cast the filtered reads as a Pandas dataframe
    detailed_results_df = pd.DataFrame.from_dict(detailed_results, orient='index')

    # Reduce this dataframe to only unique miRNA/UMI:barcode combos
    # Filter out any reads that are too far from the cropping site
    dedup_filt_df = detailed_results_df.drop_duplicates(subset=['CB_UB','miRNA']) 
    dedup_filt_df = dedup_filt_df[np.abs(dedup_filt_df['Dist_to_DROSHA']) <= Drosha_tolerance]
    
    #Generate a simple cell x gene matrix:
    count_df = dedup_filt_df.loc[:,['miRNA','CB_UB','CB']].pivot_table(
        index='miRNA', 
        columns='CB', 
        values='CB_UB', 
        aggfunc=lambda x: len(set(x))
    ).fillna(0).astype('int')

    # Make a per-miRNA tally
    counts_per_miRNA = dedup_filt_df.groupby('miRNA').agg({'UB':len}).sort_values(by='UB', ascending=False)

    print("Finished count_miRNAs: ", time.asctime())
    return detailed_results_df, count_df, counts_per_miRNA

##DEPRECATED??
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

    sample_name, bam_read_length = get_bam_features(args.BAMFILE)
    print(f'Detected sample name as: {sample_name}')

    print('Reading in gff3 file ...')
    miRNA_anno_df = read_mirbase_gff3(args.REFERENCE_FILE)
    coord_dict = make_Drosha_coord_dict(miRNA_anno_df)

    print('Detecting if BAM uses \'chr\' prefix ...')
    coord_dict = fix_chromnames(coord_dict, args.BAMFILE)

    ('Counting miRNA processing products ...')
    results, miRNA_dist_dict, counts_per_miRNA = count_miRNAs(
        args.BAMFILE, 
        coord_dict, 
        sample_name, 
        bam_read_length)

    detailed_results_df, count_df, counts_per_miRNA2 = count_miRNAs_simple(
        BAM = args.BAMFILE,
        coord_dict = coord_dict,
        sample_name = sample_name,
        bam_read_length = bam_read_length,
        Drosha_tolerance = 4, 
        write_miRNA_bam = True)
    
    # Write counts per miRNA summary
    #####

    # Distance to Drosha site pileup
    with open(os.path.join(outdir, "distance_to_DROSHA_dict.json"), "w") as outfile:
        json.dump(miRNA_dist_dict, outfile)
    #####
    print('Done with part 1!')

    raw_feature_bc_matrix = sc.read_10x_h5(args.raw)
    raw_feature_bc_matrix.var_names_make_unique()
    raw_with_miRNAs = miRNA_to_featureMatrix(results, raw_feature_bc_matrix)
    raw_with_miRNAs.var['sample_name'] = sample_name


    outfile = os.path.join(outdir, 'raw_feature_matrix_with_miRNAs.h5ad')
    raw_with_miRNAs.write(outfile)

    #v2 replacement:
    if args.results_table:
        print('Writing csv matrix of miRNA-associated reads:')
        count_df.to_csv(
            os.path.join(outdir, 'miRNA_count_matrix.csv.gz'), 
            compression='gzip', 
            index=True)

    
    if args.results_table:
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
    print(sys.argv[1:])
    main(sys.argv[1:])
