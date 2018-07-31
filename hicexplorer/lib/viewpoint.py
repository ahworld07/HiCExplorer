import numpy as np

import logging
log = logging.getLogger(__name__)
import copy

import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Viewpoint():

    def __init__(self, pHiCMatrix=None):
        self.hicMatrix = pHiCMatrix

    def readReferencePointFile(self, pBedFile):
        '''
        This function reads a text file which contains reference points. Reference points are 

        Reference points need to be tab seperated
        and contain the chromosome, the reference start point and/or the reference end point, and the gene:

        chr start gene

        or 

        chr start end gene

        Per line one reference point.

        Returns a list of tuples with (chr, start, end), if only a start index is given, end = start and a list of genes. 

        '''
        viewpoints = []
        gene_list = []
        with open(pBedFile, 'r') as file:
            for line in file.readlines():
                _line = line.strip().split('\t')
                if len(line) == 0:
                    continue
                if len(_line) == 3:
                    chrom, start, end = _line[0], _line[1], _line[1]
                    log.debug('_line: {}'.format(_line))
                    gene_list.append(_line[2])

                else:
                    chrom, start, end = _line[:3]
                    gene_list.append(_line[3])

                viewpoints.append((chrom, start, end))

        return viewpoints, gene_list

    def readInteractionFile(self, pBedFile, pRange=None):
        '''
        Reads an interaction file produced by chicViewpoint. Contains header information, these lines
        start with '#'. 
        Interactions files contain:
        Chromosome Viewpoint, Start, End, Gene, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, z-score based on relative interactions.

        This function returns:
        - header as  a string
        - interaction data in relation to relative position as a dict e.g. {-1000:0.1, -1500:0.2}  
        - rbz-score in relation to relative position as a dict (same format as interaction data)
        - interaction_file_data: the raw line in relation to the relative position. Needed for additional output file.
        '''
        # use header info to store reference point, and based matrix
        interaction_data = {}
        z_score = {}
        interaction_file_data = {}
        with open(pBedFile) as fh:
            header = fh.readline()
            for line in fh.readlines():
                # Addition header information for end users
                if line.strip().startswith('#'):
                    continue

                _line = line.strip().split('\t')
                # relative postion and relative interactions
                if int(_line[-4]) < pRange[0] or int(_line[-4]) > pRange[1]:
                    continue
                interaction_data[int(_line[-4])] = float(_line[-3])
                z_score[int(_line[-4])] = float(_line[-2])
                interaction_file_data[int(_line[-4])] = _line
        return header, interaction_data, z_score, interaction_file_data

    def readBackgroundDataFile(self, pBedFile):
        '''
        Reads a background data file, containing per line a tab delimited content:
        Relative position to viewpoint, relative interaction count to total number of interactions of all viewpoints over all samples, SEM value of this data point.
        '''
        distance = {}
        with open(pBedFile) as fh:
            for line in fh.readlines():
                _line = line.split('\t')
                distance[int(_line[0])] = [float(_line[1]), float(_line[2])]

        return distance

    def writeInteractionFile(self, pBedFile, pData, pHeader, pZscoreData):
        '''
        Writes an interaction file for one viewpoint and one sample as a tab delimited file with one interaction per line.
        Header contains information about the interaction:
        Chromosome Viewpoint, Start, End, Chromosome Interation, Start, End, Relative position (to viewpoint start / end),
        Relative number of interactions, z-score based on relative interactions, raw interaction data
        '''
        with open(pBedFile + '.bed', 'w') as fh:
            fh.write('#{}\n'.format(pHeader))
            for j, interaction in enumerate(pData):
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.12f}\t{:.12f}\t{:.12f}\n".
                         format(interaction[0], interaction[1], interaction[2],
                                interaction[3], interaction[4], interaction[5],
                                interaction[6], interaction[7], interaction[8], pZscoreData[j], interaction[9]))
        return

    def computeViewpoint(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        This function computes a viewpoint for a given sample and a given pReferencePoint within the  
        range of pRegion_start and pRegion_end. 

        All interactions with the reference point of one relative distance to it are summed up,
        if the reference point is larger than one bin of the Hi-C matrix, it is considered as one bin and the values are summed together.
        '''
        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(pReferencePoint)
        # log.debug('view_point_start {}'.format(view_point_start))
        # log.debug('view_point_end {}'.format(view_point_end))

        view_point_range = self.getViewpointRangeAsMatrixIndices(pChromViewpoint, pRegion_start, pRegion_end)
        view_point_range = list(view_point_range)
        view_point_range[1] += 1
        # log.debug('key view_point_range {} {}'.format(pRegion_start, pRegion_end))

        elements_of_viewpoint = (view_point_range[1] - view_point_range[0])
        # log.debug('view_point_range {}'.format(view_point_range))
        # log.debug('elements_of_viewpoint {}'.format(elements_of_viewpoint))
        data_list = np.zeros(elements_of_viewpoint)
        _view_point_start = view_point_start
        # TODO: check border handling! --> view_point_range[1] + 1 issue
        while _view_point_start <= view_point_end:
            chrom, start, end, _ = self.hicMatrix.getBinPos(_view_point_start)
            for j, idx in zip(range(elements_of_viewpoint), range(view_point_range[0], view_point_range[1], 1)):
                data_list[j] += self.hicMatrix.matrix[_view_point_start, idx]

            _view_point_start += 1

        elements_of_viewpoint = elements_of_viewpoint - (view_point_end - view_point_start)
        data_list_new = np.zeros(elements_of_viewpoint)
        # log.debug('elements_of_viewpoint {}'.format(elements_of_viewpoint))

        index_before_viewpoint = view_point_start - view_point_range[0]

        # elements before the viewpoint
        data_list_new[0:index_before_viewpoint] = data_list[0:index_before_viewpoint]

        # summation because the viewpoint can not be only one bin but can contain multiple
        data_list_new[index_before_viewpoint] = np.sum(data_list[index_before_viewpoint: index_before_viewpoint + view_point_end - view_point_start + 1])

        # elements after the viewpoint
        data_list_new[index_before_viewpoint + 1:] = data_list[index_before_viewpoint + view_point_end - view_point_start + 1:]
        return data_list_new

    def createInteractionFileData(self, pReferencePoint, pChromViewpoint, pRegion_start, pRegion_end, pInteractionData, pInteractionDataRaw, pGene):
        '''
        Creates out of internal information a list of tuples which can be written to an interaction file.
        Tuple contains:
        Chromosome viewpoint, start, end, chromosome interaction, start, end, relative_position, interaction data
        '''
        view_point_start, view_point_end = self.getReferencePointAsMatrixIndices(pReferencePoint)
        view_point_range = self.getViewpointRangeAsMatrixIndices(pChromViewpoint, pRegion_start, pRegion_end)
        view_point_range = list(view_point_range)
        view_point_range[1] += 1
        elements_of_viewpoint = view_point_range[1] - view_point_range[0] - (view_point_end - view_point_start) + 1

        interactions_list = []
        chrom, start, _, _ = self.hicMatrix.getBinPos(view_point_start)
        _, _, end, _ = self.hicMatrix.getBinPos(view_point_end)
        interaction_positions = list(range(view_point_range[0], view_point_start, 1))
        interaction_positions.extend([view_point_start])
        interaction_positions.extend(list(range(view_point_end + 1, view_point_range[1], 1)))
        relative_position = -1
        for j, idx in zip(range(elements_of_viewpoint), interaction_positions):

            chrom_second, start_second, end_second, _ = self.hicMatrix.getBinPos(idx)
            if relative_position < 0:
                relative_position = int(start_second) - int(start)
            else:
                relative_position = int(end_second) - int(end)

            interactions_list.append((chrom, start, end, pGene, chrom_second, start_second, end_second, relative_position, float(pInteractionData[j]), float(pInteractionDataRaw[j])))

        return interactions_list

    def getViewpointRangeAsMatrixIndices(self, pChromViewpoint, pRegion_start, pRegion_end):
        '''
        Returns the matrix indices of a chromosome and a specific position.
        '''

        _range = self.hicMatrix.getRegionBinRange(pChromViewpoint, pRegion_start, pRegion_end)
        # log.debug('_range {}'.format(_range))
        # log.debug('invert range[0], {} '.format(self.hicMatrix.getBinPos(_range[0])))
        # log.debug('invert range[1], {} '.format(self.hicMatrix.getBinPos(_range[1])))

        return _range

    def getReferencePointAsMatrixIndices(self, pReferencePoint):
        '''
        Returns the correct matrix indices of a given reference point.
        '''
        if len(pReferencePoint) == 2:
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[1]))
        elif len(pReferencePoint) == 3:
            view_point_start, view_point_end = self.hicMatrix.getRegionBinRange(pReferencePoint[0], int(pReferencePoint[1]), int(pReferencePoint[2]))
        else:
            log.error("No valid reference point given. {}".format(pReferencePoint))
            exit(1)
        return view_point_start, view_point_end

    def smoothInteractionValues(self, pData, pWindowSize):
        '''
        Adds -pWindowsSize/2 and +pWindowsSize/2 around pData[i] and averages pData[i] by pWindowSize to 
        smooth the interaction values.
        '''
        window_size = np.int(np.floor(pWindowSize / 2))
        window_size_upstream = window_size
        if pWindowSize % 2 == 0:
            window_size_upstream -= 1

        average_contacts = np.zeros(len(pData))

        # add upstream and downstream, handle regular case
        for i in range(window_size_upstream, len(pData) - window_size):
            start = i - window_size_upstream
            end = i + window_size + 1
            average_contacts[i] = np.mean(pData[start:end])

        # handle border conditions
        for i in range(window_size):
            start = i - window_size_upstream
            if start < 0:
                start = 0
            end = i + window_size + 1

            average_contacts[i] = np.mean(pData[start:end])
            average_contacts[-(i + 1)] = np.mean(pData[-end:])
        return average_contacts

    def computeRelativeValues(self, pData):
        '''
        Computes the relative values of pData by adding all data points together and divding all data points by the result.
        pData[i] = pData[i] / sum(pData)
        '''
        sumValue = np.sum(pData)
        pData /= sumValue
        return pData

    def calculateViewpointRange(self, pViewpoint, pRange):
        '''
        This function computes the correct start and end position of a viewpoint given the viewpoint and the range.
        '''
        # log.debug('self.hicMatrix.getChrBinRange(pViewpoint[0]) {}'.format(self.hicMatrix.getChrBinRange(pViewpoint[0])))
        max_length = self.hicMatrix.getBinPos(self.hicMatrix.getChrBinRange(pViewpoint[0])[1] - 1)[2]

        region_start = int(pViewpoint[1]) - pRange[0]
        if region_start < 0:
            region_start = 0

        region_end = int(pViewpoint[2]) + pRange[1]
        if region_end > max_length:
            # -1 is important, otherwise self.hicMatrix.getRegionBinRange will crash
            region_end = max_length - 1
        # log.debug('pViewpoint {}'.format(pViewpoint))
        # log.debug('pRange {}'.format(pRange))
        # log.debug('region_start {}'.format(region_start))
        # log.debug('region_end {}'.format(region_end))
        # log.debug('max_length {}'.format(max_length))

        return region_start, region_end

    def getDataForPlotting(self, pInteractionFile, pRange, pBackgroundModel):
        header, interaction_data, z_score_data, _interaction_file_data_raw = self.readInteractionFile(pInteractionFile)
        matrix_name, viewpoint, upstream_range, downstream_range, gene = header.split('\t')

        data = []
        z_score = []
        interaction_file_data_raw = {}
        if pRange:

            interaction_data_keys = copy.deepcopy(list(interaction_data.keys()))
            for key in interaction_data_keys:
                if key >= -pRange[0] and key <= pRange[1]:
                    continue
                interaction_data.pop(key, None)
            if pBackgroundModel:
                background_data_keys_sorted = sorted(pBackgroundModel)

                for key in background_data_keys_sorted:
                    if key >= -pRange[0] and key <= pRange[1]:
                        continue
                    pBackgroundModel.pop(key, None)
                background_data_keys_sorted = sorted(pBackgroundModel)

        if pBackgroundModel:
            viewpoint_index = background_data_keys_sorted.index(0)

            data_background = []
            data_background_mean = []

            for key in background_data_keys_sorted:
                if key in interaction_data:
                    data.append(interaction_data[key])

                    if key in z_score_data:
                        z_score.append(z_score_data[key])

                    data_background.append(pBackgroundModel[key][0])
                    data_background_mean.append(pBackgroundModel[key][1])

                    if key in _interaction_file_data_raw:
                        line_data_raw = _interaction_file_data_raw[key]
                        line_data_raw.append(pBackgroundModel[key][0])
                        line_data_raw.append(pBackgroundModel[key][1])

                        interaction_file_data_raw[key] = line_data_raw

        else:
            data = []
            interaction_key = sorted(interaction_data)
            for key in interaction_key:
                data.append(interaction_data[key])
            viewpoint_index = interaction_key.index(0)

        log.debug('rbz-score {}'.format(z_score))
        return header, data, data_background, data_background_mean, z_score, interaction_file_data_raw, viewpoint_index

    def plotViewpoint(self, pAxis, pData, pColor, pLabelName):
        data_plot_label = pAxis.plot(range(len(pData)), pData, '--' + pColor, alpha=0.9, label=pLabelName)

        return data_plot_label

    def plotBackgroundModel(self, pAxis, pBackgroundData, pBackgroundDataMean):
        pBackgroundData = np.array(pBackgroundData)
        pBackgroundDataMean = np.array(pBackgroundDataMean)
        data_plot_label = pAxis.plot(range(len(pBackgroundData)), pBackgroundData, '--r', alpha=0.5, label='background model')
        pAxis.fill_between(range(len(pBackgroundData)), pBackgroundData + pBackgroundDataMean, pBackgroundData - pBackgroundDataMean, facecolor='red', alpha=0.3)
        return data_plot_label

    def plotZscore(self, pAxis, pAxisLabel, pZscoreData, pLabelText, pCmap, pFigure):

        _z_score = np.empty([2, len(pZscoreData)])
        _z_score[:, :] = pZscoreData
        pAxis.xaxis.set_visible(False)
        pAxis.yaxis.set_visible(False)
        img = pAxis.contourf(_z_score, cmap=pCmap)
        divider = make_axes_locatable(pAxisLabel)
        cax = divider.append_axes("left", size="20%", pad=0.09)
        colorbar = pFigure.colorbar(img, cax=cax, ticks=[min(pZscoreData), max(pZscoreData)])

        colorbar.ax.set_ylabel('rbz-score', size=6)

        pAxisLabel.text(0.45, 0, pLabelText, size=7)
        pAxisLabel.xaxis.set_visible(False)
        pAxisLabel.yaxis.set_visible(False)
        pAxisLabel.set_frame_on(False)

    def writePlotData(self, pInteractionFileDataRaw, pFileName, pBackgroundModel):
        interaction_file_data_raw_sorted = sorted(pInteractionFileDataRaw)
        with open(pFileName + '.bed', 'w') as output_file:
            output_file.write('#ChrViewpoint\tStart\tEnd\tGene\tChrInteraction\tStart\tEnd\tRelative position\tRelative Interactions\trbz-score\tRaw')

            if pBackgroundModel:
                output_file.write('\tbackground model\tbackground model SEM\n')
            else:
                output_file.write('\n')
            for key in interaction_file_data_raw_sorted:

                _array = pInteractionFileDataRaw[key]

                output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(_array[0], _array[1], _array[2],
                                                                                      _array[3], _array[4], _array[5],
                                                                                      _array[6], _array[7], _array[8], _array[9], _array[10]))
                if pBackgroundModel:
                    output_file.write('\t{}\t{}\n'.format(_array[11], _array[12]))
                else:
                    output_file.write('\n')

    def interactionBackgroundData(self, pBackground, pRange):

        background_model = []
        background_model_sem = []

        for key in pBackground:
            if key >= -pRange[0] and key <= pRange[1]:
                background_model.append(pBackground[key][0])
                background_model_sem.append(pBackground[key][1])
                log.debug('key background {}'.format(key))
        return np.array(background_model), np.array(background_model_sem)

    def rbz_score(self, pRelativeInteractions, pBackgroundModel, pBackgroundModelSEM):
        _rbz_score = np.empty(len(pRelativeInteractions))
        if len(pRelativeInteractions) != len(pBackgroundModel) or \
                len(pRelativeInteractions) != len(pBackgroundModelSEM):
            sys.exit('Computing of rbz-score failed, data is having different size. ' +
                     '\nrelative interactions {} background model {} background model SEM {}'.format(len(pRelativeInteractions),
                                                                                                     len(pBackgroundModel), len(pBackgroundModelSEM)))
            return
        _rbz_score = pRelativeInteractions - pBackgroundModel
        _rbz_score /= pBackgroundModelSEM

        return _rbz_score
