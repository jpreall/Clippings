"""
basic IO  helper module for Clippings
"""

import collections
from collections import defaultdict
import time
import sys
import os
import shutil
import gzip
import glob
import pkg_resources
import scipy.sparse as sp
from scipy import io
import numpy as np
import h5sparse
import h5py
import pandas as pd

from scipy.sparse import csr_matrix

# From Scanpy - Fabien Theis

def _decode(item):
    try:
        item = item.decode()
    except (UnicodeDecodeError, AttributeError):
        pass

    if isinstance(item, np.ndarray):
        if item.dtype.kind == 'S':
            for dtype in ['str','int','float']:
                try:
                    item = item.astype(dtype)
                except (ValueError, UnicodeDecodeError, AttributeError):
                    pass
    return item

def _collect_datasets(group: h5py.Group, decode=True):
    dsets = {}
    
    def dig(group):
        for key, item in group.items():
            if isinstance(item, h5py.Dataset):
                if decode:
                    dsets[key] = _decode(item[:])
                else:
                    dsets[key] = item[:]
            else:
                dig(item)

    dig(group)
    
    return dsets
            
def _collect_attrs(h5attrs, decode=True):
    attrs = {}
    for key, item in h5attrs.items():
        
        if decode:
            item = _decode(item)
        
        attrs[key] = item
    
    return attrs

def read_v3_10x_h5(filename, decode=True):
    """
    Read hdf5 file from Cell Ranger v3 or later versions. 
    Updated to add h5 attributes to adata.uns 
    Adds obs key with library_id for aggregated files.
    jpreall 2022
    """
    with h5py.File(filename, 'r') as f:
        try:
            d = f['matrix']
            a = f.attrs
            dsets = _collect_datasets(d, decode=decode)
            attrs = _collect_attrs(a, decode=decode)

            M, N = dsets['shape']
            data = dsets['data']
            if dsets['data'].dtype == np.dtype('int32'):
                data = dsets['data'].view('float32')
                data[:] = dsets['data']
            matrix = sparse.csr_matrix(
                (data, dsets['indices'], dsets['indptr']),
                shape=(N, M),
                )

            adata = anndata.AnnData(
                matrix,
                obs=dict(obs_names=dsets['barcodes'].astype(str)),
                var=dict(var_names=dsets['name'].astype(str),
                            gene_ids=dsets['id'].astype(str),
                            feature_types=dsets['feature_type'].astype(str),
                            genome=dsets['genome'].astype(str),
                        ),
                uns=attrs,
                    )
            # Import library_ids into its own key in adata.obs. 
            aggr_dict = {str(n+1):key for n, key in enumerate(attrs['library_ids'])}
            adata.obs['Library_ID'] = adata.obs_names.str.split('-').str[-1].map(aggr_dict)
        except KeyError:
            raise Exception('File is missing one or more required datasets. Was this created with Cellranger < 3.0 ?')

def _read_universal_10x_h5(filename, decode=True):
    """
    Read hdf5 file from Cell Ranger v3 or later versions. Updated jpreall 2022
    """
    with h5py.File(filename, 'r') as f:
        try:
            a = f.attrs
            attrs = _collect_attrs(a, decode=decode)
            dsets = _collect_datasets(f, decode=decode)

            M, N = dsets['shape']
            data = dsets['data']

            if dsets['data'].dtype == np.dtype('int32'):
                data = dsets['data'].view('float32')
                data[:] = dsets['data']
            matrix = sparse.csr_matrix(
                (data, dsets['indices'], dsets['indptr']),
                shape=(N, M),
                )

            for key in ['gene_names','name']:
                if key in dsets.keys():
                    var_names = dsets[key]

            for key in ['genes','id']:
                if key in dsets.keys():
                    gene_ids = dsets[key]

            if 'feature_type' in dsets.keys():
                feature_types=dsets['feature_type'].astype(str),
            else:
                feature_types='Gene Expression'

            if 'genome' in dsets.keys():
                genome = dsets['genome'].astype('str')
            else:
                genome = '_'.join(list(f.keys()))

            adata = anndata.AnnData(
                matrix,
                obs=dict(obs_names=dsets['barcodes'].astype(str)),
                var = dict(var_names=var_names,
                           gene_ids=gene_ids,
                           feature_types=feature_types,
                           genome=genome),
                uns=attrs,
                    )
            # Import library_ids into its own key in adata.obs. 
            aggr_dict = {str(n+1):key for n, key in enumerate(attrs['library_ids'])}
            adata.obs['Library_ID'] = adata.obs_names.str.split('-').str[-1].map(aggr_dict)
            return adata
        
        except KeyError:
            raise Exception('File is missing one or more required datasets.')

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

def write_cellranger_h5_NEW(
    data: csr_matrix,
    template,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    library_id: str,
    genome: str = 'Unknown Genome',
    ):

    def encode_array(iterable):
        return np.array(iterable).astype('S')

    FEATURE_IDS = encode_array(var.iloc[:,0]) ## ENCODE
    FEATURE_NAMES = encode_array(var.iloc[:,1]) ## ENCODE
    FEATURE_TYPE = encode_array(var.iloc[:,2]) ##ENCODE

    ## DO TONS OF STUFF ##

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
        features.create_dataset('name', data=FEATURE_NAMES)

        f.attrs['chemistry_description'] = CHEMISTRY.encode()
        f.attrs['filetype'] = 'matrix'
        f.attrs['library_ids'] = LIBRARY_ID.encode()
        f.attrs['original_gem_groups'] = ORIG_GEM_GROUPS
        f.attrs['version'] = 2
        f.close()

def write_cellranger_h5(
    #adata,
    output_folder,
    file_prefix,
    #layer='counts',
    #aggr_key=None,
    genome=None,
    #feature_types='feature_types',
    #force_library_id=None,
    ):

    import h5py
    from scipy.sparse import csc_matrix
    from scipy.sparse import csr_matrix
    import h5sparse
    
    if layer in adata.layers.keys():
        MATRIX=csr_matrix(adata.layers[layer])
    elif layer == 'X':
        MATRIX=csr_matrix(adata.X)
        #MATRIX=csc_matrix(adata.X)
    else:
        print("Error: specified matrix layer can't be found in your data")
        
    ## identify, if it exists, the column in .var containing ensembl gene ids:
    dfstr = adata.var.select_dtypes(include=['O'])
    ens_columns = []
    for col in dfstr.columns:
        mostly_ENS_names = len([x for x in dfstr[col] if re.match('^ens', x, re.IGNORECASE)]) / len(dfstr) > 0.9
        if mostly_ENS_names:
            ens_columns += [col]
    if len(ens_columns) > 1:
        print('Warning, more than one potential ENSEMBL name column detected.  Defaulting to',ens_columns[0])
        ensembl_column = ens_columns[0]
        gene_ids = adata.var[ensembl_column]
         
    elif len(ens_columns) == 1:
        print('Using column',ens_columns[0],'as ENSEMBL ids')
        ensembl_column = ens_columns[0]
        gene_ids = adata.var[ensembl_column]
         
    else:
        print('Warning, no ENSEMBL gene id column detected.  Duplicating contents of adata.var.index')
        gene_ids = np.array(adata.var.index)

        
    ## if not already present, add 'feature_types' column to .var
    if 'feature_types' not in adata.var.columns:
        adata.var['feature_types'] = feature_types
    FEATURE_TYPE = np.array(adata.var['feature_types']).astype('S')
        
    BCS = np.array(adata.obs.index).astype('S')
    
    FEATURES = np.array(adata.var.index).astype('S')
    
    ## Detect genome  #######
    if genome:
        GENOME=genome
    else:
        if 'genome' not in adata.uns.keys():
            detect_genome(adata)

        GENOME=adata.uns['genome']
            #else:
            #    print('No genome specified, writing attribute as unspecified_genome')
            #    GENOME='unspecified_genome'
                
    ################################
        
    FEATURE_IDS = np.array(gene_ids).astype('S') 
       
    all_tag_keys = np.array([b'genome'])
    
    ## Generate list of  Library_IDs based on barcode suffixes:
    if force_library_id == None:
        sampnames = np.array(['Sample_' + suffix for suffix in adata.obs.index.str.replace('[ATGC]+-','').unique()]).astype('S')      
        LIBRARY_IDS = sampnames
        print('Auto generating Library_ID...')
        
    elif type(force_library_id) == str:
        LIBRARY_IDS = np.array([force_library_id]).astype('S')
        print('using Library_ID = ', LIBRARY_IDS)
        
    elif type(force_library_id) == list:
        LIBRARY_IDS = np.array(force_library_id).astype('S') 
        print('using Library_ID = ', LIBRARY_IDS)  
    else:
        print('Warning, force_library_id must be a str or list.')
        
    
    ORIG_GEM_GROUPS = np.array([1])

    SHAPE=adata.T.shape
    #SHAPE=adata.shape
    
    outfile=os.path.join(output_folder,file_prefix +'.h5')
    
    if aggr_key:
        ## Check if aggr_key is a categorical dtype:
        categorical_obs_cols = adata.obs.select_dtypes(include=['category']).keys()
        assert aggr_key in categorical_obs_cols, 'aggr_key must be a key corresponding to a categorical observation in .obs dataframe. Eg. \'Sample\''
        if aggr_key in categorical_obs_cols:
            obs = pd.DataFrame(adata.obs.loc[:,aggr_key])
        else:
            print('Try assigning',aggr_key,'to a categorical data type.')
            
    ## Test if barcodes are 10X Genomics compatible:
        BARCODES_IN_10X_FORMAT = len([x for x in obs.index.tolist() if re.match('[ATGCN]{16}-[0-9]+$', x, re.IGNORECASE)]) == len(adata.obs)
        
        ## Export an aggr.csv file to serve as a key for sample names:
        if BARCODES_IN_10X_FORMAT:
            obs = pd.DataFrame(adata.obs.loc[:,aggr_key])

        else:
            print("Warning, barcodes not in 10X format.  Attempting to convert.")

        ## Convert barcodes to a cellranger compatible format:
        sampnames = obs[aggr_key].unique().tolist()
        samp_num_dict = {sample:sampnames.index(sample) + 1 for sample in sampnames}
        
        ## Extract the 16 base cell barcode
        current_barcodes = obs.index.tolist()
        cell_barcode_seq = [re.search('[ATGCN]{16}', name, re.IGNORECASE).group(0).upper() for name in current_barcodes] 
        obs['cell_barcode'] = cell_barcode_seq
        obs['aggr_index'] =  obs[aggr_key].copy()
        obs = obs.replace({'aggr_index':samp_num_dict})
        obs['aggr_index'] = obs['cell_barcode'] + '-' + obs['aggr_index'].astype('str')
        obs.set_index('aggr_index', inplace=True)
        obs = pd.DataFrame(obs.loc[:,aggr_key])
        BCS = np.array(obs.index).astype('S')

        
        def write_aggr_key(obs):
            aggr_filename = os.path.join(output_folder,file_prefix+'_aggr.csv')
            print('Writing aggregation key file:',aggr_filename,'...')
            barcode_suffix = obs.index.str.replace('[ATGC]+-','')
            sampnames = obs[aggr_key]
            sample_dict = sorted(dict(zip(barcode_suffix,sampnames)).items())
            aggr_csv = pd.DataFrame(sample_dict)
            aggr_csv.columns = ['Barcode_Suffix','library_id']
            aggr_csv.to_csv(aggr_filename, index=None)
            print('Done')
        

        
        write_aggr_key(obs)
    print('Writing to',outfile)
    
    with h5sparse.File(outfile, 'w') as h5f:
        h5f.create_dataset('matrix/', data=MATRIX, compression="gzip")
        h5f.close()
    
    with h5py.File(outfile, 'r+') as f:
        f.create_dataset('matrix/barcodes', data=BCS)
        f.create_dataset('matrix/shape', (2,),dtype='int32', data=SHAPE)
        features = f.create_group('matrix/features')
        features.create_dataset('_all_tag_keys', (1,),'S6', data=all_tag_keys)  
        features.create_dataset('feature_type', data=np.array([b'Gene Expression'] * SHAPE[0]), compression="gzip")
        features.create_dataset('genome', data=np.array([GENOME.encode()] * SHAPE[0]), compression="gzip")
        features.create_dataset('id', data=FEATURE_IDS, compression="gzip")
        features.create_dataset('name', data=FEATURES, compression="gzip")
        
        f.attrs['chemistry_description'] = b'Single Cell 3\' v3'
        f.attrs['filetype'] = 'matrix'
        f.attrs['library_ids'] = LIBRARY_IDS
        f.attrs['original_gem_groups'] = ORIG_GEM_GROUPS
        f.attrs['version'] = 2
        f.close()

def read_table(file, is_gzip=False, delimiter=None):
    """
    Avoids dependency on Pandas, but this is probably dumb and overcautious
    """
    import gzip
    
    file_types = {'tsv':'\t','csv':','}
    for suffix in file_types:
        if suffix in file:
            delimiter = file_types[suffix]
            
    if file.endswith('.gz'):
        is_gzip = True

    result = []

    if is_gzip:
        f = gzip.open(file, 'r')
    else:
        f = open(file, 'r')

    try:
        for line in f:
            if isinstance(line, bytes):
                line = line.decode("utf-8")
            line = line.rstrip()
            if line:
                result.append(line.split(delimiter))
    finally:
        f.close()

    return result 

def layer_deg_and_miRNA_matrices(
    outdir, # Global output directory
    MATRIX, # Degradation matrix
    barcodes, # Full list of raw barcodes, from degradation counting step
    ):
    """
    Combine miRNA and degradation matrices to produce a unified degradation matrix with 
    0-counts for all miRNA gene entries

    MATRIX: scipy.sparse.csr_matrix of degradation counts
    barcodes: sorted list of droplet barcodes
    """
    from scipy import sparse

    # Step 1: Read in miRNA features file,
    miRNA_features_file = os.path.join(outdir,'miRNA','feature_bc_matrix_with_miRNAs','features.tsv.gz')
    all_features = read_table(miRNA_features_file)

    # Step 2: create matrix of zeros of shape: (#miRNAs x #Barcodes) (or vice versa????)
    N_NEWROWS = len(all_features) - MATRIX.shape[1]

    MIRNA_COUNTS_FILL = sparse.csr_matrix(
        np.zeros(
            [len(barcodes),N_NEWROWS]
            )
        ).astype('int32')

    # Step 4: stack deg and miRNA matrices
    COMBINED_DEG_MATRIX = sparse.hstack(
        [MATRIX,MIRNA_COUNTS_FILL],
        format='csr',
        ).astype('int32')

    return COMBINED_DEG_MATRIX, all_features
