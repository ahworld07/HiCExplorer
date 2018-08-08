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
    # parser.add_argument('--ctcfFile', '-ctcf',
    #                     help='A / B compartment file',
    #                     required=True)

    parser.add_argument('--outFileName', '-o',
                        help='File name to save the resulting matrix',
                        required=True)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    parser.add_argument('--region', '-r',
                        nargs=2,
                        help='Limit the region to x y')
    return parser


def readDomainFile(pDomainFile):

    tads = []
    with open(pDomainFile, 'r') as file:
        for i, line in enumerate(file):
            chrom, start, end, _id, score, _, _, _, _ = line.strip().split('\t')

            if i == 0:
                value_left = -1000
            else:
                value_left = tads[i - 1].valueRight

            clusterNode = ClusterNode(chrom, int(start), int(end), pValueRight=float(score), pValueLeft=value_left, pId=i)
            tads.append(clusterNode)
    return tads


def readPcaFile(pPcaFile):
    pca = []
    with open(pPcaFile, 'r') as file:
        for line in file:
            chrom, start, end, value = line.strip().split('\t')
            pca.append([chrom, start, end, value])

    return pca


def print_to_bash(pClusterNode):
    # chromosome = pChromosome
    #     self.start = pStart
    #     self.end = pEnd
    #     self.parent = pParent
    #     self.childLeft = pChildLeft
    #     self.childRight = pChildRight
    #     self.valueRight = pValueRight
    #     self.valueLeft = pValueLeft

    if pClusterNode.childLeft is None:
        return

    log.debug('id {} parent_id {} chromosome {} start {} end {} valueLeft {} valueRight {}'.format(pClusterNode.id, pClusterNode.parentId, pClusterNode.chromosome, pClusterNode.start,
                                                                                                   pClusterNode.end, pClusterNode.valueLeft,
                                                                                                   pClusterNode.valueRight))

    print_to_bash(pClusterNode.childLeft)

    print_to_bash(pClusterNode.childRight)

    # log.debug('id {} parent_id {} chromosome {} start {} end {} valueLeft {} valueRight {}'.format(pClusterNode.id, pClusterNode.parentId, pClusterNode.chromosome, pClusterNode.start,
    #                                                                                 pClusterNode.end, pClusterNode.valueLeft,
    #                                                                                 pClusterNode.valueRight))


def main(args=None):
    args = parse_arguments().parse_args(args)

    pca_data = readPcaFile(args.compartmentFile)
    tads_data = readDomainFile(args.domainsFile)
    hic_matrix = hm.hiCMatrix()

    pca_tree, _ = hic_matrix.intervalListToIntervalTree(pca_data, True)
    # domain_tree, _ = hic_matrix.intervalListToIntervalTree(domain_data, True)

    # log.debug('pca_data {}'.format(pca_tree))
    # log.debug('domain_data {}'.format(domain_tree))

    compartment_split = []
    preferences = []
    last_value = 0
    log.debug('tads_data {}'.format(tads_data[0].start))
    log.debug('pca_tree {}'.format(pca_tree['chr1'][135045000]))
    #
    _candidate_cluster = []
    node_ids = len(tads_data)
    for i in range(len(tads_data)):
        if i < len(tads_data) - 1:
            # TODO case handling if start pos is not matching with bin border of pca

            pca_value1 = list(pca_tree['chr1'][tads_data[i].start])  # [0].data
            # log.debug('pca-value1 {}'.format(pca_value1))
            if len(pca_value1) == 0:
                _candidate_cluster.append(tads_data[i])

                log.debug("pca1 zero size {}".format(tads_data[i].start))
                continue
            else:
                pca_value1 = pca_value1[0].data

            pca_value2 = list(pca_tree['chr1'][tads_data[i + 1].start])  # [0].data
            # log.debug('pca-value2 {}'.format(pca_value2))

            if len(pca_value2) == 0:
                continue
                log.debug("pca2 zero size {}".format(tads_data[i + 1].start))
            else:
                pca_value2 = pca_value2[0].data

            # log.debug('pca_valiie1 {}'.format(pca_value2))
            if np.sign(pca_value1) == np.sign(pca_value2):
                _candidate_cluster.append(tads_data[i])
            else:
                # log.debug('_candidate_cluster {}'.format(_candidate_cluster))
                compartment_split.append(_candidate_cluster)
                _candidate_cluster = []

    clustering_finished = False

    while not clustering_finished:

        merge_ids = []
        new_parent_nodes = []
        for split in compartment_split:
            preferences = []
            _merged_ids = []
            _new_parent_nodes = []
            if len(split) == 1:
                merge_ids.append(_merged_ids)
                new_parent_nodes.append(_new_parent_nodes)
                continue
            for i in range(len(split)):
                # 1 is right element, 0 is left element
                _preference = 1

                if i > 0 and i < len(split) - 1:
                    difference_left = abs(abs(split[i].valueLeft) - abs(split[i - 1].valueRight))
                    difference_right = abs(abs(split[i].valueRight) - abs(split[i + 1].valueLeft))
                    if difference_left < difference_right:
                        _preference = 0
                elif i == len(split) - 1:
                    _preference = 0

                preferences.append(_preference)
            log.debug('preferences {}'.format(preferences))
            # _merged_ids = []
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
        log.debug('First clustering phase done!')
        log.debug('mege_ids {}'.format(merge_ids))
        log.debug('new_parent_nodes {}'.format(new_parent_nodes))

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

    for cluster in compartment_split:
        log.debug('cluster parent node {}'.format(cluster))
        if len(cluster) != 0:
            print_to_bash(cluster[0])
    # log.debug("compartment_split {}".format(compartment_split))
    # for chromosome in pca_tree:
    #     for i, element in enumerate(pca_tree[chromosome]):

    #         same_compartment_segment = []
    #         if np.sign(element[2]) == np.sign(last_value):
    #             same_compartment_segment.append(domain_tree[chromosome][i])
    #             log.debug('same compartment')
    #         else:
    #             compartment_split.append(same_compartment_segment)
    #             same_compartment_segment = []
    #             same_compartment_segment.append(domain_tree[chromosome][i])

    #         last_value = element[2]

    # log.debug("compartment_split {}".format(compartment_split))
    # for chromosome in domain_tree:
    #     for i in range(len(domain_tree[chromosome])):
    #         # 1 is right element, 0 is left element
    #         _preference = 1
    #         if i > 0 and i < len(domain_tree[chromosome]) - 1:
    #             difference_left = abs(abs(domain_tree[chromosome][i][2]) - abs(domain_tree[chromosome][i - 1][2]))
    #             difference_right = (abs(domain_tree[chromosome][i][2]) - abs(domain_tree[chromosome][i + 1][2]))
    #             if difference_left < difference_right:
    #                 _preference = 0
    #         preferences.append(_preference)
    #     delete_index = []
    #     j = 0
    #     while j < range(len(preferences)):
    #         if j < len(preferences) - 1:
    #             if preferences[j] == preferences[j+1]:
    #                 clusters[j] = [[j, j+1]]
    #                 domain_tree[chromosome][j][1] = domain_tree[chromosome][j+1][1]
    #                 domain_treechromosome][j][2] = domain_tree[chromosome][j+1][2]
    #                 delete_index.append(j+1)

    #     # del dct[key]
