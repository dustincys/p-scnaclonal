#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: utils.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-04-02 20:45:33
#       History:
# =============================================================================
'''


import numpy as np
from scipy.special import gammaln
from scipy.misc import comb

from constants import EMPIRI_AAF, EMPIRI_BAF


def get_loga(data):
    return np.log(data.tumor_reads_num + 1) - np.log(data.normal_reads_num + 1)


def log_poisson_likelihood(k, Lambda):
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)


def get_cn_allele_config(max_copynumber):
    cn_allele_config = {}
    for cn in range(0, max_copynumber + 1):
        allele_config = {}
        for M_num in range(0, (cn + 2)/2):
            P_num = cn - M_num
            if P_num == 0 and M_num == 0:
                mu_T = EMPIRI_BAF / (EMPIRI_AAF + EMPIRI_BAF)
                pi_T = 'NULL'
            elif P_num == M_num:
                mu_T = 0.5
                pi_T = 'P'*P_num + 'M'*M_num
            else:
                mu_T = (M_num * 1.0) / cn
                pi_T = 'P'*P_num + 'M'*M_num + '/' + 'P'*M_num + 'M'*P_num
            allele_config[pi_T] = mu_T
        cn_allele_config[cn] = allele_config
    return cn_allele_config


def get_mu_E_joint(mu_N, mu_G, c_N, c_H, phi):
    return ((1.0 - phi)*c_N*mu_N + phi*c_H*mu_G)/((1.0 - phi)*c_N + phi*c_H)


def log_binomial_likelihood(k, n, mu):
    column_shape = (k.size, 1)
    row_shape = (1, mu.size)
    cb = comb(n, k)
    cb = cb.reshape(column_shape) * np.ones((row_shape))
    k = k.reshape(column_shape)
    n = n.reshape(column_shape)
    mu = mu.reshape(row_shape)
    return np.log(cb) + k * np.log(mu) + (n - k) * np.log(1 - mu)


def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh
