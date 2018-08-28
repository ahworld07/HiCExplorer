###
from __future__ import division
import argparse
from hicexplorer import HiCMatrix as hm
from hicexplorer._version import __version__
import numpy as np
from collections import OrderedDict
from .lib import ClusterNode
import logging
log = logging.getLogger(__name__)
from copy import deepcopy
from .utilities import expected_interactions_in_distance
from scipy.sparse import csr_matrix, linalg
from hicexplorer.utilities import toString
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros
from scipy.stats import normaltest, f_oneway, mannwhitneyu
from scipy.signal import savgol_filter

def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes hierarchical TADs given the TAD boundary file, the TAD score, A/B compartments ')

    # define the arguments
    parser.add_argument('--matrix', '-m',
                        help='interaction matrix',
                        required=True)
    # parser.add_argument('--domainsFile', '-df',
    #                     help='domains file of TADs',
    #                     required=True)

    # parser.add_argument('--compartmentFile', '-cf',
    #                     help='A / B compartment file',
    #                     required=True)
    # parser.add_argument('--ctcfFile', '-ctcf',
    #                     help='CTCF data file',
    #                     required=True)
    # parser.add_argument('--ctcfThreshold', '-ct',
    #                     help='CTCF peak threshold. If a CTCF peak is higher than is value it determines a border, the neighbor is not '
    #                     'considered as a potential nesting structure.',
    #                     required=True,
    #                     type=float)
    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix',
                        required=True,
                        nargs='+')
    parser.add_argument('--chromosomes',
                           help='List of chromosomes to be included in the '
                           'nested TAD computation.',
                           default=None,
                           nargs='+')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--region', '-r',
                        nargs=2,
                        help='Limit the region to x y')
    return parser



def writeDomainsFile(pList, pName):
    with open(pName, 'w') as file:
        for element in pList:
            file.write('{}\t{}\t{}\tx\t-1\t.\t{}\t{}\t1,2,3\n'.format(
                element[0], element[1],
                element[2], element[1],
                element[2]
            ))


def _sum_per_distance(pSum_per_distance, pData, pDistances, pN, pDistanceThreshold):
    # list_of_zero = []
    for i in range(pN):
        if not np.isnan(pData[i]):
            if pDistances[i] > pDistanceThreshold:
                continue
            pSum_per_distance[pDistances[i]] += pData[i]
            # if pDistances[i] == 0:
            #     list_of_zero.append(pData[i])
    return pSum_per_distance


def compute_zscore_matrix(pMatrix):

    instances, features = pMatrix.nonzero()

    pMatrix.data = pMatrix.data.astype(float)
    data = pMatrix.data.tolist()
    distances = np.absolute(instances - features)
    sigma_2 = np.zeros(pMatrix.shape[0])
    sum_per_distance = np.zeros(pMatrix.shape[0])
    distance_threshold = 50
    sum_per_distance = _sum_per_distance(sum_per_distance, data, distances, len(instances),pDistanceThreshold=distance_threshold)
    # compute mean
    max_contacts = np.array(range(pMatrix.shape[0], 0, -1))
    mean = sum_per_distance / max_contacts

    # compute sigma squared
    for i in range(len(instances)):
        if distances[i] > distance_threshold:
            continue
        if np.isnan(data[i]):
            sigma_2[distances[i]] += np.square(mean[distances[i]])
        else:
            sigma_2[distances[i]] += np.square(data[i] - mean[distances[i]])

    sigma_2 /= max_contacts
    sigma = np.sqrt(sigma_2)

    for i in range(len(instances)):
        if distances[i] > distance_threshold:
            data[i] = 0
            continue
        if np.isnan(data[i]):
            data[i] = (0 - mean[distances[i]]) / sigma[distances[i]]
        else:
                
            data[i] = (pMatrix.data[i] - mean[distances[i]]) / sigma[distances[i]]
            # if distances[i] < 10:
            #     log.debug('data[i] {}'.format(data[i]))
    return csr_matrix((data, (instances, features)), shape=(pMatrix.shape[0], pMatrix.shape[1]))


def expected_area(pExpected_interactions, pStart, pEnd):
    distance = pEnd - pStart
    expected_area = [pExpected_interactions[0:distance+1]] * distance
    return np.array(expected_area)
    # for i in range(distance):

def main(args=None):
    args = parse_arguments().parse_args(args)
    hic_matrix = hm.hiCMatrix(pMatrixFile=args.matrix)

    chromosome_list = []
    if args.chromosomes is not None:
        chromosome_list = args.chromosomes
    else:
        chromosome_list = hic_matrix.getChrNames()
    # pca_data = readPcaFile(args.compartmentFile, chromosome_list)
    # tads_data = readDomainFile(args.domainsFile, chromosome_list)
    # # ctcf_data = readCTCFFile(args.ctcfFile, chromosome_list)

    # # pca_tree, _ = hic_matrix.intervalListToIntervalTree(pca_data, True)

    # # compute expected values per genomic distance
    # length_chromosome = hic_matrix.matrix.shape[0]
    # chromosome_count = len(chromosome_list)
    # for chromosome in chromosome_list:

    for chromosome in chromosome_list:
        chr_range = hic_matrix.getChrBinRange(chromosome)
        log.debug('chr_range {}'.format(chr_range))
        log.debug('chr_range[0] {}'.format(chr_range[0]))
        log.debug('type chr_range {}'.format(type(chr_range[0])))

        slice_start = int(chr_range[0])

        numberOfEigenvectors = 5
        vecs_list = []
        chrom_list = []
        start_list = []
        end_list = []
        log.debug('Computing z-score matrix')
        z_score_matrix = compute_zscore_matrix(hic_matrix.matrix)
        # mask = z_score_matrix < 0
        # z_score_matrix[mask] = 0
        # z_score_matrix *= 1000
        log.debug('Computing z-score matrix...DONE')
        log.debug('Computing eigenvectors')
        window_size_pca = 50

        while slice_start < chr_range[1]:
            submatrix = z_score_matrix[slice_start:slice_start+window_size_pca, slice_start:slice_start+window_size_pca]
            corrmatrix = np.cov(submatrix.todense())
            corrmatrix = convertNansToZeros(csr_matrix(corrmatrix)).todense()
            corrmatrix = convertInfsToZeros(csr_matrix(corrmatrix)).todense()
            eigenvalues, eigenvectors = linalg.eigs(corrmatrix, k=numberOfEigenvectors)

            chrom, start, end, _ = zip(*hic_matrix.cut_intervals[slice_start:slice_start+window_size_pca])
            vecs_list += eigenvectors[:, :numberOfEigenvectors].tolist()

            chrom_list += chrom
            start_list += start
            end_list += end
            slice_start += window_size_pca
        # if args.geneTrack:
        #     vecs_list = correlateEigenvectorWithGeneTrack(ma, vecs_list, args.geneTrack)

        vecs_list = np.array(vecs_list).T
        # #sliding window to smooth values
        # window_size = 3
        threshold_window = 0.1
        for eigenvector in vecs_list:
            i = 0
            while i < (len(eigenvector)):
                if np.absolute(eigenvector[i]) < threshold_window:
                   eigenvector[i] = 0
                i += 1

        # vecs_list = np.array(vecs_list).real.T
        # vecs_list = np.array(vecs_list).T
        for i in range(len(vecs_list)):
            vecs_list[i] = savgol_filter(vecs_list[i], 5, 4)

        vecs_list = np.array(vecs_list).real.T

        for idx, outfile in enumerate(args.outFileName):
            assert(len(vecs_list) == len(chrom_list))

            with open(outfile, 'w') as fh:
                for i, value in enumerate(vecs_list):
                    if len(value) == numberOfEigenvectors:
                        if isinstance(value[idx], np.complex):
                            value[idx] = value[idx].real
                        fh.write("{}\t{}\t{}\t{:.12f}\n".format(toString(chrom_list[i]), start_list[i], end_list[i], value[idx]))

        #### detect tads: for every sign flip a tad border is reached
        vecs_list = np.array(vecs_list).T
        # threshold_eigenvector_value = []
        # for eigenvector in vecs_list:
        #     threshold_eigenvector_value.append(np.mean(np.absolute(eigenvector)) * 0.2)
        # threshold_sign = 0.1
        

        tad_list = []
        for j, eigenvector in enumerate(vecs_list):
            i = 1
            chrom_start_tad = chrom_list[0]
            start_tad = start_list[0]
            end_tad = None
            while i < len(eigenvector) - 1:
                if np.sign(eigenvector[i]) == np.sign(eigenvector[i+1]):
                    i += 1
                    continue
                if eigenvector[i] == 0:
                    i += 1
                    continue
                # if np.absolute(eigenvector[i]) < threshold_eigenvector_value[j] or np.absolute(eigenvector[i+1]) < threshold_eigenvector_value[j] :
                #     i += 1
                #     continue
                end_tad = start_list[i+1]
                tad_list.append([chrom_start_tad, start_tad, end_tad])
                chrom_start_tad = chrom_list[i+1]
                start_tad = start_list[i+1]
                i += 1


        # filter to small tads
        min_tad_size = 10

        tad_list_filtered = []
        for tad in tad_list:
            start, end = hic_matrix.getRegionBinRange(tad[0], tad[1], tad[2])
            if (end - start) < min_tad_size:
                continue
            tad_list_filtered.append(tad)
        # expected_interactions = expected_interactions_in_distance(pLength_chromosome=(chr_range[1] - chr_range[0]), pChromosome_count=21, pSubmatrix=submatrix)

        # for 
        # log.debug('vecs_list {}'.format(vecs_list))
        # z_score_theshold = 
        # p_value_threshold = 0.00001
        # tad_list_filtered = []
        # for tad in tad_list:

            # start, end = hic_matrix.getRegionBinRange(tad[0], tad[1], tad[2])
            # tad_area = hic_matrix.matrix[start:end, start:end].toarray().flatten()
            # expected = expected_area(expected_interactions, start, end).flatten()
            # # log.debug('start {} end {}'.format(start, end))
            # # log.debug('len tad_area {}  len expected {}'.format(len(tad_area),  len(expected)))
            # # log.debug('tad_area[0] {} {}'.format(tad_area[0],len(tad_area[0])))

            # # log.debug('tad_area {}'.format(tad_area))
            # test_result = mannwhitneyu(tad_area, expected)
            # pvalue = test_result[1]
            # if pvalue < p_value_threshold:
            #     tad_list_filtered.append(tad)
            # else:
            #     log.debug('rejected {} p-value {}'.format(tad, pvalue))
        writeDomainsFile(tad_list_filtered, 'tad_domains_ZscorePCA.bed')
        

    # z_score_matrix[np.logical_not(mask)] += 1
    # z_score_matrix *= z_score_matrix
    # z_score_matrix.to
    # new_z_score_matrix = []
    # for i in range(50):
        
    #     new_z_score_matrix.append(np.append(z_score_matrix.diagonal(k=i), [0]*i))
    # new_z_score_matrix = csr_matrix(new_z_score_matrix)

    hic_matrix_save = hm.hiCMatrix()

    hic_matrix_save.setMatrix(z_score_matrix, hic_matrix.cut_intervals)

    hic_matrix_save.save('nested_tads_zscore.cool')
    # corrmatrix = np.cov(new_z_score_matrix)
    # evals, eigs = linalg.eig(new_z_score_matrix.todense())


    # ctcf_tr
    # log.debug('tads_data {}'.format(tads_data))