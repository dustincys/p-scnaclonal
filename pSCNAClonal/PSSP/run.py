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

import shutil

from pSCNAClonal.PSSP.densities import FragmentSampledDensity
from pSCNAClonal.PSSP.sampler import DirichletProcessSampler
from pSCNAClonal.PSSP.trace import TraceDB

import pickle as pkl


def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data = load_data_pSCNAClonal(args.inputFilePath)

    # trace the record
    mutation_names = [data.segments[i].name for i in range(len(data.segments))]
    trace_db = TraceDB(args.out_dir, mutation_names)
    # 此处应该传递几何baseline
    # self._baseline = baseline
    # self._max_copy_number = max_copy_number
    # self._coverage = coverage

    cluster_density = FragmentSampledDensity(baseline = data.Lambda_S)

    try:
        sampler = DirichletProcessSampler(
            cluster_density,
            alpha=args.concentration,
            alpha_shape=args.concentration_prior_shape,
            alpha_rate=args.concentration_prior_rate)
    except:
        trace_db.close()
        shutil.rmtree(args.out_dir)
        raise

    sampler.sample(data.segments, trace_db, num_iters=args.num_iters)
    trace_db.close()


def load_data_pSCNAClonal(inputFilePath):
    inputFile = open(inputFilePath, 'rb')
    data = pkl.load(inputFile)
# test and debug
    data.segments = data.segments[0:100]
    return data
