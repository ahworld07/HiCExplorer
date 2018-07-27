
import pytest
from tempfile import NamedTemporaryFile, mkdtemp
import os
import shutil
from psutil import virtual_memory

from hicexplorer.utilities import genomicRegion
from hicexplorer import hicCompareMatrices

mem = virtual_memory()
memory = mem.total / 2**30

# memory in GB the test computer needs to have to run the test case
LOW_MEMORY = 2
MID_MEMORY = 7
HIGH_MEMORY = 200

REMOVE_OUTPUT = True


# matrices, zMin, zMax, --colorMap, --plotFileFormat, --plotNumbers
# method {pearson, spearman}, --log1p, --labels  --range, --outFileNameHeatmap
# --outFileNameScatter --chromosomes --threads, --help, --version
# Some definitions needed for tests
ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test_data/")

def test_help():
    args = "--help"
    hicCompareMatrices.main(args)

def test_version():
    args = "--version"
    hicCompareMatrices.main(args)