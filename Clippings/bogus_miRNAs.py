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

FILES_PATH = path.abspath(path.join(path.dirname(__file__), "../files/"))
print("This is the absolute file path: ", FILES_PATH)
print("This is pkg_resources: ", pkg_resources.resource_filename('Clippings', 'files/737K-august-2016.txt.gz'))