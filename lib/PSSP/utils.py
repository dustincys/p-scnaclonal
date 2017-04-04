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

from constants import EMPIRI_AAF, EMPIRI_BAF


def get_loga(data):
    return np.log(data.DT + 1) - np.log(data.DN + 1)


def log_poisson_likelihood(k, Lambda):
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)


def get_cn_allele_config(max_copynumber):
    cn_allele_config = {}
    for cn in range(0, max_copynumber + 1):
        allele_config = {}
        for M_num in range(0, cn + 1):
            P_num = cn - M_num
            if P_num == 0 and M_num == 0:
                mu_T = EMPIRI_BAF / (EMPIRI_AAF + EMPIRI_BAF)
                pi_T = 'NULL'
            elif P_num == M_num:
                mu_T = 0.5
                pi_T = 'P'*P_num + 'M'*M_num
            else:
                mu_T = (P_num * 1.0) / cn
                pi_T = 'P'*P_num + 'M'*M_num + '/' + 'P'*M_num + 'M'*P_num
            allele_config[pi_T] = mu_T
        cn_allele_config[cn] = allele_config
    return cn_allele_config


def get_mu_E_joint(mu_N, mu_G, c_N, c_H, phi):
    axis_1_shape = (phi.size, 1, 1)
    axis_2_shape = (1, c_H.size, 1)
    axis_3_shape = (1, 1, mu_G.size)
    phi = phi.reshape(axis_1_shape)
    c_H = c_H.reshape(axis_2_shape)
    mu_G = mu_G.reshape(axis_3_shape)
    return ((1.0 - phi)*c_N*mu_N + phi*c_H*mu_G)/((1.0 - phi)*c_N + phi*c_H)
