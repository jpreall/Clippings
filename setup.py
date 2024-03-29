import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent

VERSION = '0.1.0'
PACKAGE_NAME = 'Clippings'
AUTHOR = 'Jonathan Preall'
AUTHOR_EMAIL = 'jpreall@cshl.edu'
URL = 'https://github.com/jpreall/Clippings'

LICENSE = 'MIT'
DESCRIPTION = 'Use pysam to walk through the bam file and count every read with a ts:i tag, indicating the tso was soft clipped'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'numpy',
      'pandas',
      'scanpy',
      'matplotlib',
      'scipy',
      'seaborn',
      'pysam',
      'argparse',
      'h5py',
      'h5sparse',
      'gtfparse',
      'anndata',
      'pararead',
      'logmuse'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages(),
      entry_points={"console_scripts": ["Clippings = Clippings.run_clippings: main"]},
      )
