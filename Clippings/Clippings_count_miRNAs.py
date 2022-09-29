"""
Clippings_count_miRNAs.py

Count mapped reads harboring TSO sub-subsequence. Only works with CellRanger 4.0+ outputs.

Inputs:
    BAMFILE = BAM file produced by Cellranger v4+ with TSO reads tagged with ts:i
    REFERENCE_FILE = gff3 file with miRNA coordinates downloaded from miRBase
    outdir = path to output folder
    matrix_folder = Cellranger Count matrix folder containing gene expression data in .mtx format to merge with miRNAs
    write_bam_output (default: True) = Write candidate TSO-containing miRNA reads to BAM output file
    
    # TODO: Add ability to use miRNA GTF/coordinate files produced elsewhere
    # TODO: Change behavior of miRNA GFF3 references: Instead, specify genome version and pull from a pre-packaged file
    Usage: Clippings_count_miRNA.py BAMFILE REFERENCE_FILE --genome --outdir --raw

Author: jpreall@cshl.edu
Date: 2022-09
"""

import sys
import argparse
import os
import time
import re
import gzip
from collections import defaultdict
from typing import Union

import numpy as np
import pandas as pd
import pysam
import scipy
from scipy import io


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
    For the moment, Clippings will only look for TSO-containing reads aligning to
    the 3p-arm of annotated miRNAs. While reads piling up at the 5p-arm can sometimes
    be seen, they are more rare.
    
    Some annotated miRNAs do not have clear '3p' or '5p' arms annotated. These are
    ignored for now, but future versions will address these.
    
    
    Args:
        miRNA_anno_df (pandas DataFrame): produced by read_mirbase_gff3()
        #parent_dict (dict): Dictionary of ID to Name of only primary transcripts.

    Returns:
        coord_dict (dict): sorted by chromosome of Drosha cleavage sites of all mature miRNAs.
        Default dict has miRNA name and Drosha coordinate with '+' or '-' strand.
        
        mi_name_dict (dict): dictionary mapping miRNA gene name to "MIMAT"-style accession
    """
    
    try:
        threep = miRNA_anno_df[miRNA_anno_df['Name'].str.match('.*3p$')]
    except KeyError:
        threep = miRNA_anno_df[miRNA_anno_df['gene_name'].str.match('.*3p$')]
        
    threep = miRNA_anno_df[miRNA_anno_df['Name'].str.match('.*3p$', re.IGNORECASE)].copy()
    
    # De-duplicate any miRNAs with the same mature name, if they come from different parents
    for idx, count in threep.groupby('Name').cumcount().iteritems():
        if count != 0:
            threep.loc[idx,'Name'] = threep.loc[idx,'Name'].replace('-3p',f'.{count}-3p')

    def Drosha_site(gff_3p_row):
        COORD = gff_3p_row['start'] if \
        gff_3p_row['strand'] == '-' \
        else gff_3p_row['end']
        
        return (COORD, gff_3p_row['strand'], gff_3p_row['ID'])

    coord_dict = defaultdict(dict)
    for i, row in threep.iterrows():
        coord_dict[row.seqname][row.Name] = Drosha_site(row)
    
    mi_name_dict = dict(zip(threep['Name'],threep['ID'])) 
    
    return coord_dict, mi_name_dict
    
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

def count_miRNAs(BAM: Union[str, bytes, os.PathLike],
                        outdir: Union[str, bytes, os.PathLike],
                        coord_dict: dict, 
                        sample_name: str, 
                        bam_read_length: int,
                        Drosha_tolerance: int=4,
                        write_miRNA_bam: bool=False,
                        ):
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
        detailed_results[read_name]['miRNA'] = miRNA
        detailed_results[read_name]['CB_UB'] = read.get_tag(
            'CB') + '_' + read.get_tag('UB')
        detailed_results[read_name]['Dist_to_DROSHA'] = distance

        for tag in read.tags:
            detailed_results[read_name][tag[0]] = tag[1]
            
    def finalize_miRNA_bam(
        OUTBAM_FILENAME,
        ):
        """
        Sorts, indexes, and cleans up temporary miRNA BAM file.
        """
        print("Finalizing BAM output: ", time.asctime())    
        print("Sorting output...")
        
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
        OUTBAM_FILENAME = os.path.join(outdir,sample_name + '_sorted_miRNA_reads.bam')
        finalize_miRNA_bam(OUTBAM_FILENAME)
        
    
    ## OUTPUT 1: Cast the filtered reads as a Pandas dataframe
    detailed_results_df = pd.DataFrame.from_dict(detailed_results, orient='index')

    # Reduce this dataframe to only unique miRNA/UMI:barcode combos
    # Filter out any reads that are too far from the cropping site
    dedup_filt_df = detailed_results_df.drop_duplicates(subset=['CB_UB','miRNA']) 
    dedup_filt_df = dedup_filt_df[np.abs(dedup_filt_df['Dist_to_DROSHA']) <= Drosha_tolerance]
    
    ## OUTPUT 2: Generate a simple cell x gene matrix:
    count_df = dedup_filt_df.loc[:,['miRNA','CB_UB','CB']].pivot_table(
        index='miRNA', 
        columns='CB', 
        values='CB_UB', 
        aggfunc=lambda x: len(set(x))
    ).fillna(0).astype('int')
    
    ## OUTPUT 3: Make a table of the distance to Drosha cut site for all miRNA fragments
    ## This could be helpful in deciding which 'miRNAs' are bona-fide and which 
    ## are garbage annotations polluting miRBase, of which there are many. 
    Drosha_dist_table = dedup_filt_df.pivot_table(
        index='miRNA',
        columns=['Dist_to_DROSHA'],
        values='UB',
        aggfunc=lambda x: len(set(x))
    ).fillna(0).astype('int')

    # Extend dataframe to +/- Drosha tolerance for visual effect
    for n in list(range(-Drosha_tolerance,Drosha_tolerance+1)):
        if n not in Drosha_dist_table.columns:
            Drosha_dist_table[n] = 0

    ## OUTPUT 4: Make a per-miRNA tally
    counts_per_miRNA = dedup_filt_df.groupby('miRNA').agg({'UB':len}).sort_values(by='UB', ascending=False)

    print("Finished count_miRNAs: ", time.asctime())
    return detailed_results_df, count_df, Drosha_dist_table, counts_per_miRNA

def read_10x_mtx(matrix_path: Union[str, bytes, os.PathLike]):
    """
    Reads a Cellranger-formatted matrix path for concatentation with miRNA data
    Will accept filtered_feature_bc_matrix/ or raw_feature_bc_matrix/ folders.
    """
    VALID_10X_FILES = ['features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz']
    assert os.listdir(matrix_path) == VALID_10X_FILES, \
    f'Warning, specified Cellranger matrix folder must contain ONLY the files: {VALID_10X_FILES}'
    
    MTX = io.mmread(os.path.join(matrix_path,'matrix.mtx.gz'))
    OBS = pd.read_table(os.path.join(matrix_path,'barcodes.tsv.gz'), header=None, index_col=0)
    OBS.index.rename('Barcode', inplace=True)
    
    VAR = pd.read_table(os.path.join(matrix_path,'features.tsv.gz'), header=None, index_col=0)
    VAR.columns = ['gene_name','feature_types']
    VAR.index.rename('gene_ids', inplace=True)
    
    return MTX, OBS, VAR

def merge_miRNA_with_GEX(
    matrix_folder: Union[str, bytes, os.PathLike],
    outdir: Union[str, bytes, os.PathLike],
    count_df: pd.DataFrame ,
    mi_name_dict: dict):
    
    print(f'Reading GEX matrix from {matrix_folder}...')
    MTX, OBS, VAR = read_10x_mtx(matrix_folder)

    print(f'Merging matrices...')
    # Throw away any barcodes with miRNAs if they are not in GEX matrix
    miRNA_barcodes = set(count_df.columns)
    keep_cells = [name for name in OBS.index if name in miRNA_barcodes]
    keep_cells_indices = np.searchsorted(OBS.index, keep_cells)
    filtered_count_df = count_df.loc[:,keep_cells]
    
    # Initialize empty counts matrix to concatenate with GEX
    NEWDATA = np.zeros([len(OBS), len(count_df)], dtype='int')
    
    # Fill with miRNA counts and sparsify
    for n, row in enumerate(filtered_count_df.iterrows()):
        NEWDATA[keep_cells_indices,n] = row[1]
    sparse_mi_mtx = scipy.sparse.csr_matrix(NEWDATA.T)

    # Stack with GEX matrix
    NEWMTX = scipy.sparse.vstack([MTX,sparse_mi_mtx])

    # create feature name tsv table
    MIRVAR = pd.DataFrame(index = [mi_name_dict[mi] for mi in count_df.index])
    MIRVAR['gene_name'] = [name + '-DroshaProd' for name in count_df.index]
    MIRVAR['feature_types'] = 'Gene Expression'
    NEWVAR = pd.concat([VAR,MIRVAR])
    
    # Create outdir, if necessary
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    MTX_outdir = os.path.join(outdir,'feature_bc_matrix_with_miRNAs')

    if not os.path.exists(MTX_outdir):
        os.mkdir(MTX_outdir)
        
    # Write features.tsv.gz
    NEWVAR.to_csv(
        os.path.join(MTX_outdir,'features.tsv.gz'),
        sep = '\t',
        compression='gzip')

    print('Writing MTX..')
    # Write barcodes.tsv.gz
    OBS.to_csv(
        os.path.join(MTX_outdir,'barcodes.tsv.gz'),
        sep = '\t',
        compression='gzip',
        header=None)

    print('gzipping...')
    # Write MEX-formatted MTX file and gzip it
    MTXOUT = os.path.join(MTX_outdir,'matrix.mtx')
    io.mmwrite(
        MTXOUT,
        NEWMTX)

    MTXOUT = os.path.join(MTX_outdir,'matrix.mtx')
    with open(MTXOUT, 'rb') as FROM, gzip.open(MTXOUT + '.gz', 'wb') as TO:
        TO.writelines(FROM)
    
    # Cleanup
    os.remove(MTXOUT)
    
    file_list = [os.path.join(outdir,f) \
                 for f in ['features.tsv.gz', 'barcodes.tsv.gz','matrix.mtx.gz']]
    if all(list(map(os.path.isfile,file_list))):
        print('Successfully wrote MEX-formatted matrix with miRNAs')

def main(cmdl):
    args = _parse_cmdl(cmdl)
    outdir = args.outdir
    matrix_folder = args.matrix_folder
    #genome = args.genome

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
    coord_dict, mi_name_dict = make_Drosha_coord_dict(miRNA_anno_df)

    print('Detecting if BAM uses \'chr\' prefix ...')
    coord_dict = fix_chromnames(coord_dict, args.BAMFILE)

    #print('Counting miRNA processing products ...')
    #results, miRNA_dist_dict, counts_per_miRNA = count_miRNAs(
    #    args.BAMFILE, 
    #    coord_dict, 
    #    sample_name, 
    #    bam_read_length)

    print('Counting miRNA processing products ...')
    detailed_results_df, count_df, \
        Drosha_dist_table, counts_per_miRNA2 = count_miRNAs(
        BAM = args.BAMFILE,
        coord_dict = coord_dict,
        sample_name = sample_name,
        bam_read_length = bam_read_length,
        Drosha_tolerance = 4, 
        write_miRNA_bam = args.write_bam_output,
        outdir = outdir)
    
    ## OUTPUTS ##

    # Write detailed candidate TSO/miRNA reads in a BAM-like CSV format
    detailed_results_df.to_csv(
        os.path.join(outdir,'miRNA_read_details.csv.gz'),
        index=True,
        compression='gzip'
    )
    # Write counts per miRNA summary
    counts_per_miRNA2.to_csv(
        os.path.join(outdir,'miRNA_counts_summary.csv'),
        index=True
    )
    # Write table of TSO distances to Drosha site, for future filtering of bogus miRNAs
    Drosha_dist_table.to_csv(
        os.path.join(outdir,'TSO_distances_to_Drosha.csv'),
        index=True,
        compression='gzip')

    # Write counts matrix of miRNA processing products in all detected cells:
    print('Writing csv-formatted matrix of miRNA-associated reads:')
    count_df.to_csv(
        os.path.join(args.outdir, 'miRNA_count_matrix.csv.gz'), 
        compression='gzip', 
        index=True)

    # Merge miRNAs with existing GEX matrix
    merge_miRNA_with_GEX(
        matrix_folder = matrix_folder,
        outdir = outdir,
        count_df = count_df,
        mi_name_dict = mi_name_dict)

    #####

    # Distance to Drosha site pileup
    #with open(os.path.join(outdir, "distance_to_DROSHA_dict.json"), "w") as outfile:
    #    json.dump(miRNA_dist_dict, outfile)
    #####
    #print('Done with part 1!')

    #raw_feature_bc_matrix = sc.read_10x_h5(args.raw)
    #raw_feature_bc_matrix.var_names_make_unique()
    #raw_with_miRNAs = miRNA_to_featureMatrix(results, raw_feature_bc_matrix)
    #raw_with_miRNAs.var['sample_name'] = sample_name


    #outfile = os.path.join(outdir, 'raw_feature_matrix_with_miRNAs.h5ad')
    #raw_with_miRNAs.write(outfile)



    
    #if args.results_table:
    #    print('Writing miRNAs_result_table')
    #    results.to_csv(os.path.join(outdir, 'miRNAs_result_table.csv'), index=True)

    #print('Done with part 2!')

def _parse_cmdl(cmdl):
    """ Define and parse command-line interface. """

    parser = argparse.ArgumentParser(
        description="Clippings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('BAMFILE', help='Input bam file')
    parser.add_argument('REFERENCE_FILE', help='miRBase gff3 file')
    parser.add_argument('--outdir', dest='outdir', help='Output folder',
                        default="Clippings_outs")
    #parser.add_argument('--genome', dest='genome',
    #                    help='Genome version to record in h5 file. eg. \'hg38\' or \'mm10\'', default=None)
    #parser.add_argument('--raw', dest='raw', required=True,
    #                    help='10x Genomics raw feature bc matrix to concatenate miRNAs to')
    parser.add_argument('--matrix_folder', dest='matrix_folder', required=True,
                        help='10x Genomics matrix folder (raw or filtered) to concatenate miRNAs into')
    parser.add_argument('--write_bam_output',dest='write_bam_output',required=False, default=True,
                        help = 'Write candidate miRNA processing product reads to BAM output')

    #parser.add_argument('--results_table', dest='results_table',
    #                    help='Write out results table of reading miRNAs as csv', default=True)
    #parser = logmuse.add_logging_options(parser)
    return parser.parse_args(cmdl)

if __name__ == '__main__':
    print(sys.argv[1:])
    main(sys.argv[1:])
