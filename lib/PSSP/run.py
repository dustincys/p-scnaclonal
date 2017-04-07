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

from PSSP.densities import FragmentSampledDensity
from PSSP.sampler import DirichletProcessSampler
from PSSP.trace import TraceDB

import pickle as pkl
import numpy as np


def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data = load_data_mixclone(args.inputFilePath)
    trim_low_BAF(data)

    # trace the record
    mutation_names = [data.segments[i].name for i in range(len(data.segments))]
    trace_db = TraceDB(args.out_dir, mutation_names)
    cluster_density = FragmentSampledDensity()

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


def load_data_mixclone(inputFilePath):
    inputFile = open(inputFilePath, 'rb')
    data = pkl.load(inputFile)
    return data


def trim_low_BAF(data):
    stripes = set([seg.stripe_number for seg in data.segments])
    for sp in stripes:
        segs = filter(lambda seg: sp == seg.stripe_number, data.segments)
        trim_stripe(segs)


def trim_stripe(segments):
    coverages = []
    for seg in segments:
        coverages.extend(list(seg.BAF[:, 2] + seg.BAF[:, 3]))
    threshold = getCovThreshold(coverages)
    for seg in segments:
        d = seg.BAF[:, 2] + seg.BAF[:, 3]
        idx_rm = tuple(np.where(d < threshold)[0])
        seg.BAF = np.delete(seg.BAF, idx_rm, axis=0)


def getCovThreshold(coverages):


