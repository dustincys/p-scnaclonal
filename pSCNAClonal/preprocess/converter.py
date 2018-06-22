#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: converter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-01-31 14:07:22
#       History:
# =============================================================================
'''

import pickle as pkl
import sys
from multiprocessing import Pool

import numpy as np
import pysam

import pSCNAClonal.constants as constants
from pSCNAClonal.preprocess.data.pools.segmentPool import SegmentPool
from pSCNAClonal.preprocess.data.pools.stripePool import StripePool
from pSCNAClonal.preprocess.iofun import (PairedCountsIterator,
                                            PairedPileupIterator)
from pSCNAClonal.preprocess.mcmc import MCMCLM
from pSCNAClonal.preprocess.plotGC import GCStripePlot
from pSCNAClonal.preprocess.utils import (get_BAF_counts,
                                            normal_heterozygous_filter)

np.set_printoptions(threshold=np.inf)




class BamConverter:

    def __init__(self,
                 nBamName,
                 tBamName,
                 bedName,
                 refFaName,
                 pathPrefix="",
                 subcloneNumber=2,
                 coverage=30,
                 maxCopyNumber=6,
                 baselineThredLOH=0.3,
                 baselineThredAPM=0.01,
                 minDepth=20,
                 minBqual=10,
                 minMqual=10,
                 processNum=1,
                 bedCorrectedPath="",
                 pklPath=""):
        self._nBamName = nBamName
        self._tBamName = tBamName
        self._bedName = bedName
        self._refFaName = refFaName

        self.__pathPrefix = pathPrefix

        self.__subcloneNumber = subcloneNumber
        self.__coverage = coverage
        self.__maxCopyNumber = maxCopyNumber
        self.__baselineThredLOH = baselineThredLOH
        self.__baselineThredAPM = baselineThredAPM

        self.__minDepth = minDepth
        self.__minBqual = minBqual
        self.__minMqual = minMqual

        self.__processNum = processNum
        self.__bedCorrectedPath=bedCorrectedPath
        self.__pklPath = pklPath

        self._segPool = None

    def convert(self, readFromBed=True, method="auto", pklFlag=False):
        if not pklFlag:
            self._load_segs(readFromBed)
            # self._correct_bias(method)
            self._load_allele_counts()
            self._dump(self._segPool, ".temp.segPool")

            pklFile = open(self.__pklPath, 'wb')
            pkl.dump(self._segPool, pklFile, protocol=2)
            pklFile.close()
        else:
            pklFile = open(self.__pklPath, 'rb')
            self._segPool = pkl.load(pklFile )
            pklFile .close()

        blSegs = self._get_baseline()
        stripePool = self._generate_stripe()

        self._dump(stripePool, "stripePool.pkl")
        # self._dump_txt(stripePool, "stripePool.txt")
        self._dump(self._segPool, "segPool.pkl")

    def _dump_txt(self, stripePool, outFilePath):
        """
        out put stripePool into plain txt
        """
        stripePool.output_txt(self.__pathPrefix + outFilePath)

    def _load_allele_counts(self):
        self._get_counts(self._tBamName, self._segPool)

    def _dump(self, data, dpFileName):
        fileName = self.__pathPrefix + dpFileName
        outFile = open(fileName, 'wb')
        pkl.dump(data, outFile, protocol=2)
        outFile.close()

    def _generate_stripe(self):
        """
        generate stripe from segs
        """
        # 此处应该对每一个tag进行单独的聚类
        # 也有可能无法确定有多少个条带
        #
        # 另一种方案应该是对所有的条带进行聚类，然后从中挑选出目标条带，分解为新
        # 条带。

        yDown = constants.YDOWN
        yUp = constants.YUP
        # 此处应该近似为最大拷贝数与亚克隆数的乘积，作为手工输入也可以
        stripeNum = constants.STRIPENUML[self._segPool[-1].idx]
        noiseStripeNum = constants.NOISESTRIPENUML[self._segPool[-1].idx]

        tempSP = StripePool(self._segPool, self._segPool.baseline, yDown, yUp,
                            stripeNum, noiseStripeNum)
        tempSP.get(byTag = False)

        for idx, sp in enumerate(tempSP.stripes):
            tempSP.stripes[idx].id = idx

        return tempSP


    def _load_segs(self, readFromBed=True):
        """
        load segments for each tumor sample
        """

        # for tBamName, bedName, coverage, subcloneNumber, i in zip(self._tBamNameL,
            # self._bedNameL, self.__coverageL, self.__subcloneNumberL, range(len(self._tBamNameL))):
            # print >> sys.stdout, 'Loading segments from bam file:\n{0}\n'.format(tBamName)
            # print >> sys.stdout, 'and bed file with gc:\n{0}\n'.format(bedName)
        self._segPool = SegmentPool(self.__maxCopyNumber, coverage)
        if not readFromBed:
            nBam = pysam.Samfile(self._nBamName, 'rb')
            tBam = pysam.Samfile(tBamName, 'rb')
            self._segPool.load_seg_bam(nBam, tBam, bedName)
            nBam.close()
            tBam.close()
        else:
            self._segPool.load_seg_bed (bedName)

    def _correct_bias(self, method="auto"):
        """
        correct bias of each tumor sample
        """
        if "auto" == method:
            self._MCMC_GC_C(self._segPool, self.__subcloneNumber)
        elif "visual" == method:
            self._V_GC_C(self._segPool, len(self._segPool.segments))

    def _get_baseline(self):
        """
        get the baseline segments
        calculate baseline of each SegmentPool
        return: the baseline segments list of each SegmentPool
        """

        tempBL = segPool.get_baseline(self.__maxCopyNumber,
                                      self.__subcloneNumber,
                                      self.__baselineThredLOH,
                                      self.__baselineThredAPM,
                                      isPreprocess=True)

        return tempBL

    def _MCMC_GC_C(self, segPool, subcloneNumber):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the segPool here
        """

        mcmclm = MCMCLM(segPool, 0, subcloneNumber, self.__maxCopyNumber)
        m, c = mcmclm.run()
        print "MCMC slope = {}".format(m)

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                         np.log(seg.nReadNum + 1), segPool.segments))

        yCorrected = self._correct(x, y, m, c)

        for i in range(len(yCorrected)):
            segPool.segments[i].tReadNum = np.exp(
                yCorrected[i] +
                np.log(segPool.segments[i].nReadNum + 1)
            ) - 1
            segPool.segments[i].log_ratio = np.log(
                (yCorrected[i] + 1.0) /
                (segPool.segments[i].nReadNum + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(m, c)

    def _correct(self, x, y, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def visualize(self, segPool):
        gsp = GCStripePlot(segPool.segments, len(segPool.segments))
        print "total number: {}".format(segPool.segNum)
        gsp.plot()
        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c

    def _V_GC_C(self, segPool, sampleNumber=10000):
        gsp = GCStripePlot(segPool.segments, sampleNumber)
        print >> sys.stdout, "total number: {}".format(len(segPool.segments))
        gsp.plot()
        print >> sys.stdout, "x, y, m, c"
        print >> sys.stdout, gsp.output()

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                         np.log(seg.nReadNum + 1), segPool.segments))
        yCorrected = self._correct(x, y, m, c)

        for i in range(len(yCorrected)):
            segPool.segments[i].tReadNum = np.exp( yCorrected[i] +
                np.log(segPool.segments[i].nReadNum + 1)) - 1
            segPool.segments[i].log_ratio = np.log(
                (yCorrected[i] + 1.0) /
                (segPool.segments[i].nReadNum + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def _get_counts(self, tBamName, segPool):
        """
        get allele counts of target bam file
        save the counts into segPool
        """

        segNum = len(segPool.segments)
        processNum = self.__processNum
        print "processNum = {}".format(processNum)

        if processNum > segNum:
            processNum = segNum

        pool = Pool(processes=processNum)

        argsL = []

        for j in range(0, segNum):
            segName = segPool.segments[j].name
            chromName = segPool.segments[j].chromName
            chromIdx = segPool.segments[j].chromIdx
            start = segPool.segments[j].start
            end = segPool.segments[j].end

            argsT = (
                segName,
                chromName,
                chromIdx,
                start,
                end,
                self._nBamName,
                tBamName,
                self._refFaName,
                self.__minDepth,
                self.__minBqual,
                self.__minMqual)

            argsL.append(argsT)

        countsTL = pool.map(process_by_segment, argsL)

        for j in range(0, segNum):
            pairedCountsJ, BAFCountsJ = countsTL[j]
            segPool.segments[j].pairedCounts = pairedCountsJ
            segPool.segments[j].BAFCounts = BAFCountsJ

# ===============================================================================
#  Function
# ===============================================================================


def process_by_segment(argsT):
    segName, chromName, chromIdx, start, end, nBamName,\
        tBamName, refFaName, minDepth, minBqual,\
        minMqual = argsT

    print 'Preprocessing segment {0}...'.format(segName)
    sys.stdout.flush()

    nBam = pysam.Samfile(nBamName, 'rb')
    tBam = pysam.Samfile(tBamName, 'rb')
    refFasta = pysam.Fastafile(refFaName)

    normalPileupIter = nBam.pileup(chromName, start, end)
    tumorPileupIter = tBam.pileup(chromName, start, end)

    pairedPileupIter = PairedPileupIterator(
        normalPileupIter, tumorPileupIter, start, end)
    pairedCountsIter = PairedCountsIterator(
        pairedPileupIter,
        refFasta,
        chromName,
        chromIdx,
        minDepth,
        minBqual,
        minMqual)

    pairedCountsJ, BAFCountsJ = iterator_to_counts(pairedCountsIter)
    countsTuple_j = (pairedCountsJ, BAFCountsJ)

    nBam.close()
    tBam.close()
    refFasta.close()

    return countsTuple_j


def iterator_to_counts(pairedCountsIter):
    buff = 100000

    pairedCountsJ = np.array([[], [], [], [], [], []], dtype=int).transpose()
    BAFCountsJ = np.zeros((100, 100))
    buffCounts = []
    i = 0

    for counts in pairedCountsIter:
        buffCounts.append(counts)
        i = i + 1

        if i < buff:
            continue

        buffCounts = np.array(buffCounts)

        if buffCounts.shape[0] != 0:
            BAFCountsBuff = get_BAF_counts(buffCounts)
            BAFCountsJ += BAFCountsBuff

        buffCountsFiltered = normal_heterozygous_filter(buffCounts)

        if buffCountsFiltered.shape[0] != 0:
            pairedCountsJ = np.vstack((pairedCountsJ, buffCountsFiltered))

        buffCounts = []
        i = 0

    buffCounts = np.array(buffCounts)

    if buffCounts.shape[0] != 0:
        BAFCountsBuff = get_BAF_counts(buffCounts)
        BAFCountsJ += BAFCountsBuff

    buffCountsFiltered = normal_heterozygous_filter(buffCounts)

    if buffCountsFiltered.shape[0] != 0:
        pairedCountsJ = np.vstack((pairedCountsJ, buffCountsFiltered))

    return (pairedCountsJ, BAFCountsJ)
