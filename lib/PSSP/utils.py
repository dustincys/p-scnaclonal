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
                pi_T = ''
                mu_T = -1
            else:
                pi_T = 'P'*P_num + 'M'*M_num
                mu_T = (P_num * 1.0) / cn

            allele_config[pi_T] = mu_T
        cn_allele_config[cn] = allele_config

    return cn_allele_config
