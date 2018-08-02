###
from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes hierachical TADs given the TAD boundary file, the TAD score, A/B compartments ')

    # define the arguments
    parser.add_argument('--domainsFile', '-df',
                        help='domains file of TADs',
                        required=True)

    parser.add_argument('--compartmentFile', '-cf',
                        help='A / B compartment file',
                        required=True)
    parser.add_argument('--ctcfFile', '-ctcf',
                        help='A / B compartment file',
                        required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser

def readDomainFile(pDomainFile):

    tads = []
    with open(pDomainFile, 'r') as file:
        for line in file:
            chrom, start, end, _, score = line.strip().split('\t')
            tads.append([chrom, start, end, score])
    return tads
def readPcaFile(pPcaFile):
    pca = []
    with open(pDomainFile, 'r') as file:
        for line in file:
            chrom, start, end, value = line.strip().split('\t')
            tads.append([chrom, start, end, value])

    return pca


def main(args=None):
    args = parse_arguments().parse_args(args)

    