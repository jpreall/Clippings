"""
Combine miRNA_adata.h5 and deg_count_adata.h5 into a single .h5ad as layers
"""
import scanpy as sc
import argparse
import sys
import os
import glob

def load_from_folder(folder_path, min_counts=200, verbose=True):
    """
    Loading and concatenating multiple samples' raw miRNA with deg_count .h5ad files.
    """
    matrix_list = glob.glob(os.path.join(folder_path, '*.h5ad'))
    initial_matrix = sc.read(matrix_list[0])
    if verbose:
        print(initial_matrix.var['sample_name'][0] + ' prefilter shape: ', initial_matrix.shape)
        sc.pp.filter_cells(initial_matrix, min_counts=min_counts)
        print(initial_matrix.var['sample_name'][0] + ' post-filter shape: ', initial_matrix.shape)
    else:
        sc.pp.filter_cells(initial_matrix, min_counts=min_counts)

    non_initial_dict = {}  # dictionary holding all anndata objs to be combined, but not the initial one
    for file in matrix_list[1:]:
        non_initial_dict[file] = sc.read(file)
        if verbose:
            print()
            print(non_initial_dict[file].var['sample_name'][0] + ' prefilter shape: ', non_initial_dict[file].shape)
            sc.pp.filter_cells(non_initial_dict[file], min_counts=min_counts)
            print(non_initial_dict[file].var['sample_name'][0] + ' prefilter shape: ', non_initial_dict[file].shape)
        else:
            sc.pp.filter_cells(non_initial_dict[file], min_counts=min_counts)

    # concatenate the matrices to the initial matrix
    # set batch category equal to the sample_name of each matrix for identification
    combined = initial_matrix.concatenate(*list(non_initial_dict.values()), batch_key=None,
                               batch_categories=[initial_matrix.var['sample_name'][0]] + [
                                   list(non_initial_dict.values())[index].var['sample_name'][0] for index in
                                   range(len(list(non_initial_dict.values())))], join='outer')
    return combined


def main(args):
    outdir = args.outdir
    deg_adata = sc.read_10x_h5(args.deg_adata)
    deg_adata.var_names_make_unique()
    miRNA_adata = sc.read_h5ad(args.miRNA_adata)
    miRNA_adata.var_names_make_unique()
    print('deg_adata dimensions: ', deg_adata)
    print('miRNA_adata dimensions: ', miRNA_adata)

    print('barcodes only in miRNA and not deg: ', miRNA_adata[list(set(list(miRNA_adata.obs_names)) - set(list(deg_adata.obs_names)))])
    print('miRNA genes only in miRNA adata and not deg adata: ', miRNA_adata[:,list(set(list(miRNA_adata.var_names)) - set(list(deg_adata.var_names)))])
    non_deg_barcodes = miRNA_adata[list(set(list(miRNA_adata.obs_names)) - set(list(deg_adata.obs_names)))]
    # set all barcode values equal to zero **** this is important (need unit testing)
    non_deg_barcodes.X[non_deg_barcodes.X>0] = 0

    miRNA_for_deg_barcodes = miRNA_adata[deg_adata.obs_names,list(set(list(miRNA_adata.var_names)) - set(list(deg_adata.var_names)))].copy()
    # set all miRNA gene values equal to zero **** this is important (need unit testing)
    miRNA_for_deg_barcodes.X[miRNA_for_deg_barcodes.X>0] = 0
    print('miRNA genes: ', miRNA_for_deg_barcodes.var_names)

    # combining miRNA genes for deg barcodes with the deg adata (ie. 33694 genes + 75 miRNA genes = 33769 total genes)
    # Need unit testing
    tmp_combine = deg_adata.T.concatenate(miRNA_for_deg_barcodes.T, join='outer', index_unique=None)
    print('transpose merged dimensions: ', tmp_combine)
    deg_with_miRNAs = tmp_combine.T
    print('deg_with_miRNAs dimensions: ', deg_with_miRNAs)
    del tmp_combine

    # join barcodes from raw that are not in deg adata (end result should be same dimension as the miRNA adata in order to make layers)
    # Need unit testing
    deg_with_miRNAs_allBarcodes = deg_with_miRNAs.concatenate(non_deg_barcodes, join='outer', index_unique=None)
    print('combined and merged dimensions of deg_with_miRNAs_allBarcodes: ', deg_with_miRNAs_allBarcodes)
    print('Need to add as layer to miRNA_adata with dimensions: ', miRNA_adata)

    # match barcode and var_name order from deg_with_miRNAs_allBarcodes to miRNA_adata
    print('Need to make sure that the layers are ordered correspondingly in terms of barcodes and var_names')
    print('Start with barcode ordering')
    print('miRNA_adata barcode #50000-50003: ', miRNA_adata[50000:50003].obs.index)
    print('miRNA_adata barcode #50000-50003 values: ', miRNA_adata[50000:50003].X.data)
    print('deg_with_miRNAs_allBarcodes barcode #50000-50003: ', deg_with_miRNAs_allBarcodes[50000:50003].obs.index)
    print('deg_with_miRNAs_allBarcodes barcode #50000-50003 values: ', deg_with_miRNAs_allBarcodes[50000:50003].X.data)

    deg_with_miRNAs_allBarcodes_sorted_bar = deg_with_miRNAs_allBarcodes[miRNA_adata.obs.index]
    print('deg_with_miRNAs_allBarcodes_sorted_bar barcode #50000-50003: ', deg_with_miRNAs_allBarcodes_sorted_bar[50000:50003].obs.index)
    print('deg_with_miRNAs_allBarcodes_sorted_bar barcode #50000-50003 values: ', deg_with_miRNAs_allBarcodes_sorted_bar[50000:50003].X.data)
    print('Done with barcode matching/ordering')

    print('Move onto var_names ordering')
    print('miRNA_adata var_names #100-103: ', miRNA_adata[:,100:103].var_names)
    print('miRNA_adata var_names #100-103 values: ', miRNA_adata[:,100:103].X.data)
    print('deg_with_miRNAs_allBarcodes_sorted_bar var_names #100-103: ', deg_with_miRNAs_allBarcodes_sorted_bar[:,100:103].var_names)
    print('deg_with_miRNAs_allBarcodes_sorted_bar var_names #100-103 values: ', deg_with_miRNAs_allBarcodes_sorted_bar[:,100:103].X.data)

    deg_with_miRNAs_allBarcodes_sorted_bar_var = deg_with_miRNAs_allBarcodes_sorted_bar[:,miRNA_adata.var_names]
    print('deg_with_miRNAs_allBarcodes_sorted_bar_var var_names #100-103: ', deg_with_miRNAs_allBarcodes_sorted_bar_var[:,100:103].var_names)
    print('deg_with_miRNAs_allBarcodes_sorted_bar_var var_names #100-103 values: ', deg_with_miRNAs_allBarcodes_sorted_bar_var[:,100:103].X.data)
    print('Done with var_names matching/ordering')

    # set the layer to the matching/ordered deg_count
    miRNA_adata.layers['deg_counts'] = deg_with_miRNAs_allBarcodes_sorted_bar_var.X
    print('Finish adding deg_counts layer!')


    if os.path.isdir(outdir):
        # should throw error
        print('Output directory already exists')
    os.mkdir(outdir)
    #os.chdir(outdir)

    # should we specify .h5ad?
    outfile=os.path.join(outdir,'raw_feature_matrix_with_miRNAs_and_deg_count.h5ad')
    miRNA_adata.write(outfile)

    # testing purposes
    foo = sc.read(outfile)
    print('miRNA table dimensions: ', foo)

    #if args.results_table:
    #    print('Writing miRNAs_result_table')
    #    miRNA_adata.X.to_csv(os.path.join(outdir, 'miRNAs_result_table.csv'), index=True)
    #    print('Writing deg_counts_result_table')
    #    miRNA_adata.layers['deg_counts'].to_csv(os.path.join(outdir, 'deg_counts_result_table.csv'), index=True)
    #else:
    #    print('User did not want csv file')
    print('Done with part 2!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--deg_adata', dest='deg_adata', required=True, help='Path to degradation anndata')
    parser.add_argument('--miRNA_adata', dest='miRNA_adata', required=True, help='Path to miRNA anndata')
    parser.add_argument('--outdir', dest='outdir', help='Output folder', default='combined_deg_miRNA_adata')
    #parser.add_argument('--results_table', dest='results_table', help='Writes out layers as csvs', default=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    main(args)
