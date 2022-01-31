import gzip
import argparse
import math
import numpy as np
import argparse
import itertools
from distutils.dir_util import copy_tree
import os



parser = argparse.ArgumentParser(description='Join sequencing reads into a single UMIed sequence')
parser.add_argument('--input_dir', help='the directory to copy',required=True)
parser.add_argument('--sample', help='sample name',required=True)
parser.add_argument('--project', help='the name of this project')
parser.add_argument('--webdir', help='the base web output directory',required=True)


args = parser.parse_args()

# do we need to make the project directory
full_path = args.webdir + '/' + args.project + '/' + args.sample + '/'

copy_tree(args.input_dir, full_path)
