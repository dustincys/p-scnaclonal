'''
Created on 2012-11-21

@author: andrew
'''
from multiprocessing import Pool
from collections import namedtuple

from pydp.densities import Density, log_binomial_pdf, log_poisson_pdf
# from pydp.utils import log_sum_exp

from utils import get_loga


FragmentSampledData = namedtuple(
    'FragmentSampledData', ['CN', 'SN', 'DT', 'DN', 'BAFs', 'BASELINE'])


class FragmentSampledDensity(Density):
    def __init__(self, baseline=0, max_copy_number=7, params=None):
        Density.__init__(self, params=None)
        self.baseline = baseline
        self.max_copy_number = max_copy_number

    def log_p(self, data, params):
        # params.phi

        ll = 0.0

        copy_numbers = None

        if data.BASELINE == "True":
            copy_numbers = [2]
        elif get_loga(data) < self.baseline:
            copy_numbers = range(2, self.max_copy_number + 1)
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
        self._getRD(data, copy_number, phi)

        genotype, mu = getAlleleConfig(copy_number)

        return ll_seg, pi

    def _getRD(self, data, copy_number, phi):
        bar_c = phi * copy_number + (1.0 - phi) * 2.0
        lambda_possion = (bar_c / 2.0) * baseline * (data.DN + 1) - 1
        ll_RD = log_poisson_pdf(data.DT, lambda_possion)
        return ll_RD

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
