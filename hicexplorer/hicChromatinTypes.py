from __future__ import division

import argparse
import numpy as np
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__

import logging
log = logging.getLogger(__name__)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     add_help=False,
                                     description=('Takes two matrices as input, normalizes them and applies '
                                                  'the given operation. To normalize the matrices '
                                                  'each element is divided by the sum of the matrix.'))

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--pca', '-m',
                                help='Name of the matrices in .h5 format to use, separated by a space.',
                                metavar='matrix.h5',
                                nargs=3,
                                required=True)

    # parserRequired.add_argument('--outFileName', '-o',
    #                             help='File name to save the resulting matrix. The output is '
    #                             'also a .h5 file.',
    #                             required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--operation',
                           help='Operation to apply to the matrices.',
                           choices=['diff', 'ratio', 'log2ratio'],
                           default='log2ratio')

    parserOpt.add_argument("--help", "-h", action="help", help="show this help message and exit")

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):

    args = parse_arguments().parse_args(args)
    # if args.operation not in ['diff', 'ratio', 'log2ratio']:
