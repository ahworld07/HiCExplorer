from __future__ import division
import argparse
from past.builtins import zip
from scipy.sparse import lil_matrix

from hicexplorer.iterativeCorrection import iterativeCorrection
from hicmatrix import HiCMatrix as hm
from hicexplorer._version import __version__
from hicexplorer.utilities import toString
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from hicexplorer.utilities import check_cooler

import numpy as np
debug = 0

import logging
log = logging.getLogger(__name__)

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        conflict_handler='resolve',
        description="""


Iterative correction for a Hi-C matrix (see Imakaev et al. 2012
Nature Methods for details). For the method to work correctly, bins with
zero reads assigned to them should be removed as they can not be corrected.
Also, bins with low number of reads should be removed,
otherwise, during the correction step, the counts associated with
those bins will be amplified (usually, zero and low coverage bins
tend contain repetitive regions).  Bins with extremely high number
of reads can also be removed from the correction as they may represent
copy number variations.

To aid in the identification of bins with low and high read coverage, the
histogram of the number of reads can be plotted together with the
Median Absolute Deviation (MAD).

It is recommended to run hicCorrectMatrix as follows:

    $ hicCorrectMatrix diagnostic_plot --matrix hic_matrix.h5 -o plot_file.png

Then, after revising the plot and deciding the threshold values:

    $ hicCorrectMatrix correct --matrix hic_matrix.h5 \r
    --filterThreshold <lower threshold> <upper threshold> -o corrected_matrix

For a more in-depth review of how to determine the threshold values, please visit:
http://hicexplorer.readthedocs.io/en/latest/content/example_usage.html#correction-of-hi-c-matrix
"""
    )

    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(
        title="Options",
        dest='command',
        metavar='',
        help="""To get detailed help on each of the options: \r

    $ chicCorrectMatrix diagnostic_plot -h \r

    $ chicCorrectMatrix correct -h
    """)
    plot_mode = subparsers.add_parser(
        'diagnostic_plot',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="""Plots a histogram of the coverage per bin together with the modified
z-score based on the median absolute deviation method
(see Boris Iglewicz and David Hoaglin 1993, Volume 16: How to Detect
and Handle Outliers The ASQC Basic References in Quality Control:
Statistical Techniques, Edward F. Mykytka, Ph.D., Editor).
        """,
        usage='%(prog)s '
              '--matrix hic_matrix.h5 '
              '-o file.png')
    plot_modeRequired = plot_mode.add_argument_group('Required arguments')
    plot_modeRequired.add_argument('--matrix', '-m',
                                   help='Name of the Hi-C matrix to correct in .h5 format.',
                                   required=True)

    plot_modeRequired.add_argument('--plotName', '-o',
                                   help='File name to save the diagnostic plot.',
                                   required=True)

    plot_modeOpt = plot_mode.add_argument_group('Optional arguments')
    plot_modeOpt.add_argument('--chromosomes',
                              help='List of chromosomes to be included in the iterative '
                              'correction. The order of the given chromosomes will be then '
                              'kept for the resulting corrected matrix.',
                              default=None,
                              nargs='+')

    plot_modeOpt.add_argument('--xMax',
                              help='Max value for the x-axis in counts per bin.',
                              default=None,
                              type=float)

    plot_modeOpt.add_argument(
        '--perchr',
        help='Compute histogram per chromosome. For samples from cells with uneven number '
        'of chromosomes and/or translocations it is advisable to check the histograms '
        'per chromosome to find the most conservative `filterThreshold`.',
        action='store_true')

    plot_modeOpt.add_argument('--verbose',
                              help='Print processing status.',
                              action='store_true')

    subparsers.add_parser(
        'correct',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[correct_subparser()],
        help="""Run the iterative correction.""",
        usage='%(prog)s '
              '--matrix hic_matrix.h5 '
              '--filterThreshold -1.2 5 '
              '-out corrected_matrix.h5 \n')

    return parser


def correct_subparser():
    # define the arguments
    parser = argparse.ArgumentParser(add_help=False)

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='Name of the Hi-C matrix to correct in .h5 format.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name to save the resulting matrix. The '
                                'output is a .h5 file.',
                                required=True)

    parserRequired.add_argument('--filterThreshold', '-t',
                                help='Removes bins of low or large coverage. '
                                'Usually these bins do not contain valid Hi-C data or represent '
                                'regions that accumulate reads and thus must be discarded. '
                                'Use hicCorrectMatrix diagnostic_plot '
                                'to identify the modified z-value thresholds. A lower and upper '
                                'threshold are required separated by space, e.g. --filterThreshold '
                                '-1.5 5',
                                type=float,
                                nargs=2,
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--iterNum', '-n',
                           help='Number of iterations to compute.',
                           type=int,
                           metavar='INT',
                           default=500)

    parserOpt.add_argument('--inflationCutoff',
                           help='Value corresponding to the maximum number of times a bin '
                           'can be scaled up during the iterative correction. For example, '
                           'an inflation cutoff of 3 will filter out all bins that were '
                           'expanded 3 times or more during the iterative correction.',
                           type=float)

    parserOpt.add_argument('--transCutoff', '-transcut',
                           help='Clip high counts in the top -transcut trans '
                           'regions (i.e. between chromosomes). A usual value '
                           'is 0.05 ',
                           type=float)

    parserOpt.add_argument('--sequencedCountCutoff',
                           help='Each bin receives a value indicating the '
                           'fraction that is covered by reads. A cutoff of '
                           '0.5 will discard all those bins that have less '
                           'than half of the bin covered.',
                           default=None,
                           type=float)

    parserOpt.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the iterative '
                           'correction. The order of the given chromosomes will be then '
                           'kept for the resulting corrected matrix',
                           default=None,
                           nargs='+')

    parserOpt.add_argument('--skipDiagonal', '-s',
                           help='If set, diagonal counts are not included',
                           action='store_true')

    parserOpt.add_argument('--perchr',
                           help='Normalize each chromosome separately. This is useful for '
                           'samples from cells with uneven number of chromosomes and/or translocations.',
                           action='store_true')

    parserOpt.add_argument('--verbose',
                           help='Print processing status',
                           action='store_true')
    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))
    return parser




def main(args=None):
    args = parse_arguments().parse_args(args)
    if args.verbose:
        log.setLevel(logging.INFO)

    # args.chromosomes
    if check_cooler(args.matrix) and args.chromosomes is not None and len(args.chromosomes) == 1:
        ma = hm.hiCMatrix(args.matrix, pChrnameList=toString(args.chromosomes))
    else:
        ma = hm.hiCMatrix(args.matrix)

        if args.chromosomes:
            ma.reorderChromosomes(toString(args.chromosomes))

    # mask all zero value bins
    row_sum = np.asarray(ma.matrix.sum(axis=1)).flatten()
    log.info("Removing {} zero value bins".format(sum(row_sum == 0)))
    ma.maskBins(np.flatnonzero(row_sum == 0))
    matrix_shape = ma.matrix.shape
    ma.matrix = convertNansToZeros(ma.matrix)
    ma.matrix = convertInfsToZeros(ma.matrix)

    if 'plotName' in args:
        plot_total_contact_dist(ma, args)
        log.info("Saving diagnostic plot {}\n".format(args.plotName))
        return

    log.info("matrix contains {} data points. Sparsity {:.3f}.".format(
        len(ma.matrix.data),
        float(len(ma.matrix.data)) / (ma.matrix.shape[0] ** 2)))

    if args.skipDiagonal:
        ma.diagflat(value=0)

    outlier_regions = filter_by_zscore(ma, args.filterThreshold[0], args.filterThreshold[1], perchr=args.perchr)
    # compute and print some statistics
    pct_outlier = 100 * float(len(outlier_regions)) / ma.matrix.shape[0]
    ma.printchrtoremove(outlier_regions, label="Bins that are MAD outliers ({:.2f}%) "
                                               "out of".format(pct_outlier, ma.matrix.shape[0]),
                        restore_masked_bins=False)

    assert matrix_shape == ma.matrix.shape
    # mask filtered regions
    ma.maskBins(outlier_regions)
    total_filtered_out = set(outlier_regions)

    if args.sequencedCountCutoff and 0 < args.sequencedCountCutoff < 1:
        chrom, _, _, coverage = zip(*ma.cut_intervals)

        assert type(coverage[0]) == np.float64

        failed_bins = np.flatnonzero(
            np.array(coverage) < args.sequencedCountCutoff)

        ma.printchrtoremove(failed_bins, label="Bins with low coverage", restore_masked_bins=False)
        ma.maskBins(failed_bins)
        total_filtered_out = set(failed_bins)
        """
        ma.matrix, to_remove = fill_gaps(ma, failed_bins)
        log.warning("From {} failed bins, {} could "
                         "not be filled\n".format(len(failed_bins),
                                                  len(to_remove)))
        ma.maskBins(to_remove)
        """

    if args.transCutoff and 0 < args.transCutoff < 100:
        cutoff = float(args.transCutoff) / 100
        # a usual cutoff is 0.05
        ma.truncTrans(high=cutoff)

    pre_row_sum = np.asarray(ma.matrix.sum(axis=1)).flatten()
    correction_factors = []
    if args.perchr:
        corrected_matrix = lil_matrix(ma.matrix.shape)
        # normalize each chromosome independently
        for chrname in list(ma.interval_trees):
            chr_range = ma.getChrBinRange(chrname)
            chr_submatrix = ma.matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]]
            _matrix, _corr_factors = iterative_correction(chr_submatrix, args)
            corrected_matrix[chr_range[0]:chr_range[1], chr_range[0]:chr_range[1]] = _matrix
            correction_factors.append(_corr_factors)
        correction_factors = np.concatenate(correction_factors)

    else:
        corrected_matrix, correction_factors = iterative_correction(ma.matrix, args)

    ma.setMatrixValues(corrected_matrix)
    ma.setCorrectionFactors(correction_factors)
    log.info("Correction factors {}".format(correction_factors[:10]))
    if args.inflationCutoff and args.inflationCutoff > 0:
        after_row_sum = np.asarray(corrected_matrix.sum(axis=1)).flatten()
        # identify rows that were expanded more than args.inflationCutoff times
        to_remove = np.flatnonzero(after_row_sum / pre_row_sum >= args.inflationCutoff)
        ma.printchrtoremove(to_remove,
                            label="inflated >={} "
                            "regions".format(args.inflationCutoff), restore_masked_bins=False)
        total_filtered_out = total_filtered_out.union(to_remove)

        ma.maskBins(to_remove)

    ma.printchrtoremove(sorted(list(total_filtered_out)),
                        label="Total regions to be removed", restore_masked_bins=False)

    ma.save(args.outFileName, pApplyCorrection=False)
