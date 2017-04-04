'''
# =============================================================================
#      FileName: densities.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-04-03 13:11:45
#       History: @author: andrew
# =============================================================================
'''
from multiprocessing import Pool
from collections import namedtuple

from pydp.densities import Density, log_binomial_pdf, log_poisson_pdf
# from pydp.utils import log_sum_exp

from utils import get_loga, get_cn_allele_config

import constants


FragmentSampledData = namedtuple(
    'FragmentSampledData', ['CN', 'SN', 'DT', 'DN', 'BAFs', 'BASELINE'])


class FragmentSampledDensity(Density):

    def __init__(self, baseline=0, max_copy_number=7, params=None):
        Density.__init__(self, params=None)

        self._baseline = baseline
        self._max_copy_number = max_copy_number
        self._allele_config = get_cn_allele_config(max_copy_number)

    def log_p(self, data, params):
        # params.phi

        ll = 0.0

        copy_numbers = None

        if data.BASELINE == "True":
            copy_numbers = [2]
        elif get_loga(data) < self._baseline:
            copy_numbers = range(2, self._max_copy_number + 1)
        else:
            copy_numbers = range(0, 2 + 1)

        ll = ll + self._log_RD() + self._log_BAF()

        return ll

#       for cn_r, cn_v, mu_v, log_pi in zip(
#               data.cn_r, data.cn_v, data.mu_v, data.log_pi):
#           temp = log_pi + self._log_binomial_likelihood(data.b,
#                                                         data.d,
#                                                         params.phi,
#                                                         params.s,
#                                                         data.eps,
#                                                         mu_v,
#                                                         cn_r,
#                                                         cn_v)

#           ll.append(temp)

#       return log_sum_exp(ll)

    def _getLLSeg(self, data, copy_number, phi):
        # get the likelihood of segment and it genotype
        ll_seg = 0

        ll_rd= self._getRD(data, copy_number, phi)

        allele_types = self._allele_config[copy_number]
        ll_baf, pi = self._getBAF(data, copy_number, allele_types, phi)

        ll_seg = ll_baf + ll_rd

        return ll_seg, pi

    def _getRD(self, data, copy_number, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        bar_c = phi * copy_number + (1.0 - phi) * c_N
        lambda_possion = (bar_c / c_N) * baseline * (data.DN + 1) - 1
        ll_RD = log_poisson_pdf(data.DT, lambda_possion)
        return ll_RD

    def _getBAF(self, data, copy_number, allele_types, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        mu_N = constants.MU_N
        mu_G = allele_types.keys()
        mu_E = get_mu_E_joint(mu_N, mu_G, c_N, copy_number, phi)
        a_T_j = data.BAF[:, 2]
        b_T_j = data.BAF[:, 3]
        d_T_j = a_T_j + b_T_j
        # add prior or not?
        ll = log_binomial_likelihood_joint(b_T_j, d_T_j, mu_E)
        ll_LOH_j = np.logaddexp.reduce(ll, axis=3).sum(axis=0)
        return ll_LOH_j

    def _log_binomial_likelihood(self, b, d, phi, s, eps, mu_v, cn_r, cn_v):
        cn_N = 2

        # Compute probabilities of sampling from ref and var populations
        p_r = (1 - s) * cn_N + s * (1 - phi) * cn_r
        p_v = s * phi * cn_v

        # Normalise probabilities
        p_r = p_r / (p_r + p_v)
        p_v = 1 - p_r

        mu = p_r * eps + p_v * mu_v

        return log_binomial_pdf(b, d, mu)
