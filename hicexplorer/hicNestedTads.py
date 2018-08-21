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

def parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes hierarchical TADs given the TAD boundary file, the TAD score, A/B compartments ')

    # define the arguments
    parser.add_argument('--domainsFile', '-df',
                        help='domains file of TADs',
                        required=True)

    parser.add_argument('--compartmentFile', '-cf',
                        help='A / B compartment file',
                        required=True)
    parser.add_argument('--ctcfFile', '-ctcf',
                        help='CTCF data file',
                        required=True)
    parser.add_argument('--ctcfThreshold', '-ct',
                        help='CTCF peak threshold. If a CTCF peak is higher than is value it determines a border, the neighbor is not '
                        'considered as a potential nesting structure.',
                        required=True,
                        type=float)
    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix',
                        required=True)
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


def readDomainFile(pDomainFile, pChromosomeList):

    tads = []
    with open(pDomainFile, 'r') as file:
        for i, line in enumerate(file):
            chrom, start, end, _id, score, _, _, _, _ = line.strip().split('\t')

            if chrom in pChromosomeList:

                if i == 0:
                    value_left = -1000
                else:
                    value_left = tads[i - 1].valueRight

                clusterNode = ClusterNode(chrom, int(start), int(end), pValueRight=float(score), pValueLeft=value_left, pId=i)
                tads.append(clusterNode)
    return tads


def readPcaFile(pPcaFile, pChromosomeList):
    pca = []
    with open(pPcaFile, 'r') as file:
        for line in file:
            chrom, start, end, value = line.strip().split('\t')
            if chrom in pChromosomeList:
                pca.append([chrom, start, end, value])

    return pca

def readCTCFFile(pCtcfFile, pChromosomeList):
    ctcf = []
    with open(pCtcfFile, 'r') as file:
        for line in file:
            chrom, start, end, value = line.strip().split('\t')
            if chrom in pChromosomeList:
                ctcf.append([chrom, int(start), int(end), float(value)])

    return ctcf

def createList(pClusterNode, pList):
    if pClusterNode.childLeft is not None:
        pList.append(createList(pClusterNode.childLeft, pList))
    if pClusterNode.childRight is not None:
        pList.append(createList(pClusterNode.childRight, pList))
    return [pClusterNode.chromosome,
            pClusterNode.start,
            pClusterNode.end,
            pClusterNode.id,
            pClusterNode.valueRight,
            pClusterNode.start,
            pClusterNode.end]


def writeDomainsFile(pList, pName):
    with open(pName, 'w') as file:
        for element in pList:
            file.write('{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t1,2,3\n'.format(
                element[0], element[1],
                element[2], element[3],
                element[4], element[5],
                element[6]
            ))


def print_to_bash(pClusterNode):
    # chromosome = pChromosome
    #     self.start = pStart
    #     self.end = pEnd
    #     self.parent = pParent
    #     self.childLeft = pChildLeft
    #     self.childRight = pChildRight
    #     self.valueRight = pValueRight
    #     self.valueLeft = pValueLeft
    log.debug('id {} parent_id {} chromosome {} start {} end {} valueLeft {} valueRight {}'.format(pClusterNode.id, pClusterNode.parentId, pClusterNode.chromosome, pClusterNode.start,
                                                                                                   pClusterNode.end, pClusterNode.valueLeft,
                                                                                                   pClusterNode.valueRight))
    if pClusterNode.childLeft is not None:
        print_to_bash(pClusterNode.childLeft)

    if pClusterNode.childRight is not None:
        print_to_bash(pClusterNode.childRight)

def match_peaks_to_bin(pIntervalTree, pPeaks):


    # in which bin is the peak?
    adjusted_peaks = {}
    for peak in pPeaks:
        # log.debug('peak {}'.format(peak))
        # log.debug('{}'.format( list(pIntervalTree[peak[0]][peak[1]])[0].begin ))
        
        _data = list(pIntervalTree[peak[0]][peak[1]])
        if len(_data) == 0:
            continue
        _data = _data[0]
        if _data.begin  in adjusted_peaks:
            adjusted_peaks[_data.begin][3] += peak[3]
        else:
            adjusted_peaks[_data.begin] = [peak[0], _data.begin, _data.end, peak[3]] 
        
    return list(adjusted_peaks.values())


def main(args=None):
    args = parse_arguments().parse_args(args)
    chromosome_list = []
    if args.chromosomes is not None:
        chromosome_list = args.chromosomes
    else:
        chromosome_list = hic_matrix.getChrNames
    pca_data = readPcaFile(args.compartmentFile, chromosome_list)
    tads_data = readDomainFile(args.domainsFile, chromosome_list)
    hic_matrix = hm.hiCMatrix()
    ctcf_data = readCTCFFile(args.ctcfFile, chromosome_list)

    pca_tree, _ = hic_matrix.intervalListToIntervalTree(pca_data, True)
    # ctcf_tree = deepcopy(pca_tree)

    ctcf_tree, _ = hic_matrix.intervalListToIntervalTree(match_peaks_to_bin(pca_tree, ctcf_data), True)
    # log.debug('{}'.format(ctcf_tree))
    # domain_tree, _ = hic_matrix.intervalListToIntervalTree(domain_data, True)

    # log.debug('pca_data {}'.format(pca_tree))
    # log.debug('domain_data {}'.format(domain_tree))

    compartment_split = []
    preferences = []
    last_value = 0
    # for tad in tads_data:
    #     log.debug('tads_data {} {} '.format(tad.valueLeft, tad.valueRight))

    # log.debug('pca_tree {}'.format(pca_tree['chr1'][135045000]))
    #
    for chromosome in chromosome_list:
        _candidate_cluster = []
        node_ids = len(tads_data)
        for i in range(len(tads_data)):
            if i < len(tads_data) - 1:

                pca_value1 = list(pca_tree[chromosome][tads_data[i].start])  # [0].data
                pca_value2 = list(pca_tree[chromosome][tads_data[i + 1].start])
                if len(pca_value1) == 0 or len(pca_value2) == 0:
                    _candidate_cluster.append(tads_data[i])
                    continue
                pca_value1 = pca_value1[0].data
                pca_value2 = pca_value2[0].data

                if np.sign(pca_value1) == np.sign(pca_value2):

                    # startbin = sorted(self.interval_trees[chrname][startpos:startpos + 1])[0].data
                    # endbin = sorted(self.interval_trees[chrname][endpos:endpos + 1])[0].data
  
                    ctcf_value = list(ctcf_tree[chromosome][tads_data[i + 1].start])
                    if not len(ctcf_value) == 0:
                        if ctcf_value[0].data >= args.ctcfThreshold:
                            _candidate_cluster.append(tads_data[i])
                    continue
                
                _candidate_cluster.append(tads_data[i])
                compartment_split.append(_candidate_cluster)
                _candidate_cluster = []
            else:
                pca_value1 = list(pca_tree[chromosome][tads_data[i].start])
                pca_value2 = list(pca_tree[chromosome][tads_data[i - 1].start])
                if len(pca_value1) == 0 or len(pca_value2) == 0:
                    compartment_split[-1].append(tads_data[i])
                    continue
                pca_value1 = pca_value1[0].data
                pca_value2 = pca_value2[0].data
                if np.sign(pca_value1) == np.sign(pca_value2):
                    compartment_split[-1].append(tads_data[i])
                else:
                    compartment_split.append([tads_data[i]])

        clustering_finished = False

        while not clustering_finished:

            merge_ids = []
            new_parent_nodes = []
            for k, split in enumerate(compartment_split):
                preferences = []
                _merged_ids = []
                _new_parent_nodes = []
                if len(split) == 1:
                    merge_ids.append(_merged_ids)
                    new_parent_nodes.append(_new_parent_nodes)
                    continue
                for i, tad in enumerate(split):
                    # 1 is right element, 0 is left element
                    _preference = 0

                    if i == 0:
                        _preference = 1
                    elif i == len(split) - 1:
                        _preference = 0
                    elif tad.valueLeft < tad.valueRight:
                        # preference to right element
                        _preference = 1

                    preferences.append(_preference)

                for i in range(len(split)):
                    if i < len(split) - 1:

                        if preferences[i] == 1 and preferences[i + 1] == 0:
                            clusterNode = ClusterNode(pChromosome=split[i].chromosome,
                                                    pStart=split[i].start,
                                                    pEnd=split[i + 1].end,
                                                    pValueLeft=split[i].valueLeft,
                                                    pValueRight=split[i + 1].valueRight,
                                                    pParent=None,
                                                    pChildLeft=split[i],
                                                    pChildRight=split[i + 1],
                                                    pId=node_ids,
                                                    pChildLeftId=split[i].id,
                                                    pChildRightId=split[i + 1].id)
                            split[i].parent = clusterNode
                            split[i].parentId = node_ids
                            split[i + 1].parent = clusterNode
                            split[i + 1].parentId = node_ids

                            node_ids += 1
                            _merged_ids.append((i, i + 1))
                            _new_parent_nodes.append(clusterNode)
                merge_ids.append(_merged_ids)
                new_parent_nodes.append(_new_parent_nodes)

            for i, split in enumerate(merge_ids):
                for j, ids in enumerate(split):
                    compartment_split[i][ids[0]] = new_parent_nodes[i][j]
                    compartment_split[i][ids[1]] = None
                for ids in compartment_split[i]:
                    if ids is None:
                        compartment_split[i].remove(ids)
            counter = 0
            for split in compartment_split:
                if len(split) == 1:
                    counter += 1
            if counter == len(compartment_split):
                clustering_finished = True

            for _ids in new_parent_nodes:
                if len(_ids) == 0:
                    clustering_finished = True
                else:
                    clustering_finished = False
                    break

        cluster_list = []

        for cluster in compartment_split:
            # log.debug('cluster parent node {}'.format(cluster))
            if len(cluster) != 0:
                # print_to_bash(cluster[0])
                cluster_list.append(createList(cluster[0], cluster_list))

        # log.debug('cluster_list {}'.format(cluster_list))

        writeDomainsFile(cluster_list, args.outFileName + '_' + chromosome)
