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
from collections import namedtuple
from pydp.densities import Density, log_poisson_pdf
from utils import get_loga, get_cn_allele_config, get_mu_E_joint,\
    log_binomial_likelihood, mad_based_outlier
import constants
from random import randint
import numpy as np


SegmentResultData = namedtuple(
    'SegmentResultData', ['likelihood', 'copy_number', 'cellular_prevalence'])


class FragmentSampledDensity(Density):

    def __init__(self, baseline=0, max_copy_number=7, coverage=30, params=None):
        Density.__init__(self, params=None)

        self._baseline = baseline
        self._max_copy_number = max_copy_number
        self._coverage = coverage
        self._allele_config = get_cn_allele_config(max_copy_number)

    def log_p(self, seg, phi):
        # phi
        if seg.phi_last is not None and phi in seg.phi_last.keys():
            return seg.phi_last[phi].likelihood
        else:
            if seg.phi_last is None:
                seg.phi_last = {}

            if len(seg.phi_last) >= constants.SEG_RECORD_SIZE:
                idx = randint(0, constants.SEG_RECORD_SIZE - 1)
                seg.phi_last.pop(seg.phi_last.keys()[idx])

            ll, cn, pi = self._getSegResData(seg, phi)
            seg.phi_last[phi] = SegmentResultData(ll, cn, pi)
            return ll

    def _getSegResData(self, seg, phi):
        copy_numbers = None
        if seg.BASELINE == "True":
            copy_numbers = [2]
        elif get_loga(seg) > self._baseline:
            copy_numbers = range(2, self._max_copy_number + 1)
        else:
            copy_numbers = range(0, 2 + 1)

        ll_pi_s = [self._getLLSeg(seg, copy_number, phi) for copy_number in
                   copy_numbers]
        (ll, pi) = max(ll_pi_s, key=lambda x: x[0])
        cn = ll_pi_s.index((ll, pi))
        return ll, cn, pi

    def _getLLSeg(self, seg, copy_number, phi):
        ll_seg = 0
        ll_rd = self._getRD(seg, copy_number, phi)
        allele_types = self._allele_config[copy_number]
        self._augBAF(seg, copy_number)
        if 0 == seg.BAF.shape[0]:
            ll_baf = 0
            pi = "*"
        else:
            ll_baf, pi = self._getBAF(seg, copy_number, allele_types, phi)
        ll_seg = ll_baf + ll_rd
        return ll_seg, pi

    def _augBAF(self, seg, copy_number):
        if copy_number > 2:
            threshold = constants.BAF_THRESHOLD * self._coverage
            d_T_j = np.sum(seg.BAF[:, 2:4], axis=1)
            idx_rm = tuple(np.where(d_T_j < threshold)[0])
            seg.BAF = np.delete(seg.BAF, idx_rm, axis=0)
        else:
            pass

        if seg.BAF.shape[0] > 1:
            b_T_j = np.min(seg.BAF[:, 2:4], axis=1)
            d_T_j = np.sum(seg.BAF[:, 2:4], axis=1)
            baf = b_T_j * 1.0 / d_T_j
            outlier = mad_based_outlier(baf)
            seg.BAF = np.delete(seg.BAF, list(outlier.astype(int)), axis=0)
        else:
            pass

    def _getRD(self, seg, copy_number, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        bar_c = phi * copy_number + (1.0 - phi) * c_N
        lambda_possion = (bar_c / c_N) * self._baseline * (seg.DN + 1) - 1
        ll_RD = log_poisson_pdf(seg.DT, lambda_possion)
        return ll_RD

    def _getBAF(self, seg, copy_number, allele_types, phi):
        c_N = constants.COPY_NUMBER_NORMAL
        mu_N = constants.MU_N
        mu_G = np.array(allele_types.keys())
        mu_E = get_mu_E_joint(mu_N, mu_G, c_N, copy_number, phi)
        b_T_j = np.min(seg.BAF[:, 2:4], axis=1)
        d_T_j = np.sum(seg.BAF[:, 2:4], axis=1)
        # add prior or not?
        ll = log_binomial_likelihood(b_T_j, d_T_j, mu_E)
        ll_bafs = ll.sum(axis=0)
        idx_max = ll_bafs.argmax(axis=0)
        ll_baf = ll_bafs[idx_max]
        pi = allele_types[allele_types.keys()[idx_max]]
        return ll_baf, pi
