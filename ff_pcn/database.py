#!/usr/bin/env python2

"""
Database storing results.
"""

import os
import logging
from ff_pcn import ExistanceReasonRegular, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox, ExistanceReasonNeedFactorization


RESULT_FOLDER = os.path.abspath(os.path.join(__file__, '../../result/'))


class Database(object):

    def __init__(self, result_folder=RESULT_FOLDER):
        self.result_folder = result_folder

    def add(self, p, e, n, result):
        if isinstance(result, ExistanceReasonRegular):
            resstring = 'regular'
        elif isinstance(result, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox):
            resstring = 'L > U'
        elif isinstance(result, ExistanceReasonNeedFactorization):
            resstring = 'Factorization needed'
        with open(os.path.join(self.result_folder, 'ex_%d.txt' % n), 'a') as fp:
            fp.write('(%d, %d %d) => %s\n' % (p, e, n, resstring))

database = Database()
