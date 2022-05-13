import shutil

import gtfparse

from Clippings import deg_count_with_UMIs
from pararead import ParaReadProcessor
import pytest
import os, os.path
import pandas as pd

@pytest.fixture(scope='module')
def deg_count_parser():
    print('------------------setup-----------------')
    #db = deg_count_with_UMIs.ReadCounter(ParaReadProcessor)
    bam_file_loc = os.path.join(os.getcwd(), '500_PBMC_3p_LT_Chromium_X_possorted_mini_bam.bam')
    parser = deg_count_with_UMIs._parse_cmdl([bam_file_loc, '--TSSgtf',
                                              'GRCh38-2020-A-103genes.gtf', '--outdir', 'test_outdir', '--mtx', 'True'])
    yield parser
    print('------------------teardown------------------')


@pytest.fixture(scope='module')
def deg_count_dict_of_TSSes(deg_count_parser):
    print('------------------setup-----------------')
    TSS_dict, feature_dict = deg_count_with_UMIs.dict_of_TSSes(deg_count_parser.TSSgtf)
    yield TSS_dict, feature_dict
    print('------------------teardown------------------')


@pytest.fixture(scope='module')
def deg_count_read_counter(deg_count_parser, deg_count_dict_of_TSSes):
    print('------------------setup-----------------')
    initial_directory = os.getcwd()
    tmp_outdir = os.path.join(initial_directory, deg_count_parser.outdir)

    if not os.path.isdir(tmp_outdir):
        os.mkdir(tmp_outdir)
    os.chdir(tmp_outdir)

    # Creating directory to write out json dict
    jsonFolder = 'jsonFolder'
    jsonPath = os.path.join(os.getcwd(), jsonFolder)
    os.mkdir(jsonPath)

    counter = deg_count_with_UMIs.ReadCounter(deg_count_parser.readsfile,
                                              cores=deg_count_parser.cores,
                                              outfile=deg_count_parser.outdictionary, action="CountReads",
                                              TSSdictionary=deg_count_dict_of_TSSes[0],
                                              features=deg_count_dict_of_TSSes[1],
                                              write_degraded_bam_file=deg_count_parser.write_degraded_bam_file,
                                              include_introns=deg_count_parser.include_introns)
    counter.register_files()
    good_chromosomes = counter.run()
    counter.combine(good_chromosomes, chrom_sep="\n")
    os.chdir(initial_directory)
    yield counter
    print('------------------teardown------------------')
    shutil.rmtree(str(tmp_outdir))


#def test_file_path():



def test_parser(deg_count_parser):
    assert deg_count_parser.cores == 10
    assert deg_count_parser.readsfile == os.path.join(os.getcwd(), '500_PBMC_3p_LT_Chromium_X_possorted_mini_bam.bam')
    assert deg_count_parser.outdictionary == 'dictionary.txt'
    assert deg_count_parser.TSSgtf == 'GRCh38-2020-A-103genes.gtf'
    assert deg_count_parser.outdir == 'test_outdir'
    assert deg_count_parser.genome is None
    assert deg_count_parser.mtx == 'True'  # argparse cannot take a boolean (need string)
    assert deg_count_parser.write_degraded_bam_file is False
    assert deg_count_parser.include_introns is False


def test_dict_of_TSSes(deg_count_dict_of_TSSes):
    TSS_dict = deg_count_dict_of_TSSes[0]
    feature_dict = deg_count_dict_of_TSSes[1]
    assert len(TSS_dict) == 103  # the number of genes in the initial gtf file
    assert len(feature_dict) == 103  # the number of genes in the initial gtf file


@pytest.mark.usefixtures("deg_count_read_counter")
def test_json_chromosome_count(deg_count_read_counter, deg_count_parser):
    test_bed = pd.read_csv('500_PBMC_3p_LT_Chromium_X_possorted_mini_bed.bed', sep='\t', comment='t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
              'blockCount', 'blockSizes', 'blockStarts']
    test_bed.columns = header[:len(test_bed.columns)]
    num_bed_chroms = len(test_bed['chrom'].unique())

    jsonFolder = 'jsonFolder'
    jsonPath = os.path.join(os.getcwd(), deg_count_parser.outdir, jsonFolder)
    print('jsonPath: ', jsonPath)
    num_json_files = len([name for name in os.listdir(jsonPath) if os.path.isfile(os.path.join(jsonPath, name))])

    assert num_json_files == num_bed_chroms
    print('num_json_files: ', num_json_files)
    print('num_bed_chroms: ', num_bed_chroms)
