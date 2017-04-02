#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
# =============================================================================
#      FileName: run.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-03-28 09:12:15
# =============================================================================
'''

from __future__ import division

from math import log

import csv
import shutil

from PSSP.densities import FragmentSampledData, FragmentSampledDensity
from PSSP.sampler import DirichletProcessSampler
from PSSP.trace import TraceDB


def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data = load_data(args.inputFilePath)

    trace_db = TraceDB(args.out_dir, data.keys())
    cluster_density = FragmentSampledDensity()

    try:
        sampler = DirichletProcessSampler(
            cluster_density,
            args.tumour_content,
            alpha=args.concentration,
            alpha_shape=args.concentration_prior_shape,
            alpha_rate=args.concentration_prior_rate)
    except:
        trace_db.close()
        shutil.rmtree(args.out_dir)
        raise

    sampler.sample(data.values(), trace_db, num_iters=args.num_iters)
    trace_db.close()

def load_data(in_file):
    inputFile = open(inputFilePath, 'rb')
    data_pkl = pkl.load(inputFile)

    data = {}

    for i in range(data.pkl.seg_num):
        segName = data_pkl.segments[i]

        CN = data.segments[i].copy_number
        SN = data.segments[i].stripe_number
        DT = data.segments[i].tumor_reads_num
        DN = data.segments[i].normal_reads_num
        BAFs = data.segments[i].BAF_counts
        BASELINE = data.segments[i].baseline_label

        data[segName] = FragmentSampledData(CN,
                                            SN,
                                            DT,
                                            DN,
                                            BAFs,
                                            BASELINE)

    return data



