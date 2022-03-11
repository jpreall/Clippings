import sys
import os.path as path
import pkg_resources
import pandas as pd
import scanpy as sc


def find_families(adata):
    """
    Label families for given AnnData in .var section and create dataframe for miRNA
    families and their associated counts.

    A bit longer description.

    Args:
        adata (AnnData): AnnData object with microRNA genes added on.
        Normally read in from .h5ad created from Clippings_count_miRNAs.py

    Returns:
        pandas DataFrame: DataFrame containing list of miRNA families and the
        associated counts for each family.

    Raises:
        Exception: None

    """
    FILES_PATH = path.abspath(path.join(path.dirname(__file__), "../files/"))
    print("This is the absolute file path: ", FILES_PATH)
    print("This is pkg_resources: ", pkg_resources.resource_filename(
        'Clippings', 'files/multi_mi_familes.csv'))
    families = pd.read_csv(pkg_resources.resource_filename(
        'Clippings', 'files/multi_mi_familes.csv'))

    # filter min counts for cells and genes to ensure 'n_counts' label is found
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_counts=1)

    # remove duplicate groups of miRNAs (usually a +/- strand difference, but sometimes they are in multiple genes)
    # in the future, there is likely a better way to do this, but currently the family names are not important
    families_dedup = families.drop_duplicates(subset=['miRNAs']).copy()

    # turn the miRNAs in each family into a list
    families_dedup['miRNAs'] = families_dedup['miRNAs'].str.split(';')

    families_dedup['combined'] = families_dedup['Chromosome'].astype(str) + ':' + families_dedup['start'].astype(
        str) + '-' + families_dedup['end'].astype(str)

    # since the number of miRNAs in explode before and after dropping duplicates is slightly different, this means
    # there are some individual miRNAs in multiple families
    fam_dedup_explode = families_dedup.explode('miRNAs').drop_duplicates(subset=['miRNAs'])
    micro_fam_dict = dict(zip(fam_dedup_explode['miRNAs'], fam_dedup_explode['Gene_name']))
    adata.var['family'] = adata.var['gene_ids'].map(micro_fam_dict)
    adata.var['family'] = adata.var['family'].fillna(adata.var['gene_ids'])
    family_count = adata[:, adata.var_names.str.contains('_DroshaProd')].var.groupby(['family']).\
        n_counts.sum().reset_index().sort_values('n_counts', ascending=False)
    return family_count
