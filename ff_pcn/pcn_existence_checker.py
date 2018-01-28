#!/usr/bin/env sage

"""
Module holding classes finding PCNs on computational basis.
"""

__author__ = "Stefan Hackenberg"

try:
    import ff_pcn
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(__file__, '../../')))
import logging
import multiprocessing
from sage.all import Integer, euler_gamma, log, PolynomialRing, divisors, factor, primes, ZZ
from ff_pcn import ExistanceReasonRegular, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox, ExistanceReasonNeedFactorization
from ff_pcn.basic_number_theory import is_regular
from ff_pcn.finite_field_extension import FiniteFieldExtension
from ff_pcn.database import database


def check_p_n(pn):
    p, n = pn
    for e in xrange(1,n):
        q = p**e
        if q > n:
            break
        checker = PCNExistenceChecker(p, e, q, n)
        res = checker.check_existance()
        logging.getLogger(__name__).info('check_until_n of (%d, %d, %d) => %s', p, e, n, res)
        database.add(p, e, n, res[1])
        del res
        del checker


def check_n_multiprocessing(n):
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(check_p_n, ((p,n) for p in primes(n)))


def check_n(n):
    for p in primes(n):
        check_p_n((p, n))


class PCNExistenceChecker(object):
    """
    Class for checking existence of PCN element in
    E/F where F = GF(q^n) and F = GF(q) with q = p^r.
    """

    @staticmethod
    def check_range(l, m):
        """
        Checks exitance of PCNs for all FiniteField extensions of degree from l to m.
        """
        # for n in xrange(l, m):
        #     check_n(n)
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(check_n, xrange(l,m))

    @staticmethod
    def check_to(m):
        """
        Checks exitance of PCNs for all FiniteField extensions of degree to m.
        """
        results = dict()
        for n in xrange(m):
            for p in primes(n):
                for e in xrange(1,n):
                    q = p**e
                    if q > n:
                        break
                    checker = PCNExistenceChecker(p, e, q, n)
                    res = results[(p,e,n)] = checker.check_existance()
                    logging.getLogger(__name__).info('check_until_n of (%d, %d, %d) => %s', p, e, n, res)
        return results

    def __init__(self, p, e, q, n):
        self.p = p
        self.e = e
        self.n = n
        self.q = q
        self.ff_extension = FiniteFieldExtension(self.p, self.e, self.q, self.n)

    def check_existance(self):
        logging.getLogger(__name__).debug('check_existance of %s', self)
        p = self.p
        q = self.q
        n = self.n
        e = self.e

        # Existance if (p,e,n) is regular
        if is_regular(p, e, 1, n, 1):
            logging.getLogger(__name__).debug('check_existance:is regular')
            return (True, ExistanceReasonRegular(self))

        # Existance if q > n
        if (q > n):
            logging.getLogger(__name__).debug('check_existance:q > n')
            return (True, ExistanceReasonQBiggerN(self))

        # Otherwise check the more complex criterions
        result = self._check_existance_primitives_more_equal_not_normals_approx()
        if result is True:
            logging.getLogger(__name__).debug('check_existance: |P| >= |H| approx')
            return (True, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox(self))

        return (False, ExistanceReasonNeedFactorization(self))

    def _check_existance_primitives_more_equal_not_normals_approx(self):
        logging.getLogger(__name__).debug('_check_existance_primitives_more_equal_not_normals_approx')
        n = self.n
        q = self.q
        e = self.e
        qn = q*n

        # Check approximation for euler_phi:
        # euler_phi(n) >= n/( e^gamma * loglog n + 3/loglog n )
        euler_phi_approx = (qn-1)/(e**euler_gamma * log(log(qn-1)) + 3/log(log(qn-1)))
        if euler_phi_approx > self.ff_extension.upper_border_not_normals():
            return True

    def __str__(self):
        return "<PCNExistenceChecker q=%d^%d n=%d>" % (self.p, self.e, self.n)


if __name__ == '__main__':
    import sys
    from ff_pcn.datastore import datastore
    logging.basicConfig(level=logging.INFO)
    PCNExistenceChecker.check_range(int(sys.argv[1]), int(sys.argv[2]))
