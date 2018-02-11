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
from sage.all import Integer, euler_gamma, log, PolynomialRing, divisors, factor, primes, ZZ, prod, euler_phi
from ff_pcn import ExistanceReasonRegular, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox, ExistanceReasonPrimitivesMoreEqualNotNormals, ExistanceReasonNeedFactorization, ExistanceReasonFoundOne, ExistanceReasonNotExisting
from ff_pcn.basic_number_theory import is_regular, factor_with_euler_phi
from ff_pcn.finite_field_extension import FiniteFieldExtension
from ff_pcn.database import database
from ff_pcn.factorer import factorer


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


def qs_to_check(n):
    """
    Returns a list of qs which has to be checked for given n.
    """
    qs = []
    for p in primes(n):
        for e in xrange(1,n):
            q = p**e
            if q > n:
                break
            qs += [(p,e)]
    return qs


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
        #     check_n_multiprocessing(n)
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(check_n, xrange(l, m))

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

    def check_existance(self, no_explicit_search=True):
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

        self.factorization = self._factor()
        if self.factorization is None:
            return (False, ExistanceReasonNeedFactorization(self))

        result = self._check_existance_primitives_more_equal_not_normals()
        if result is True:
            logging.getLogger(__name__).debug('check_existance: |P| >= |H|')
            return (True, ExistanceReasonPrimitivesMoreEqualNotNormals(self))

        if no_explicit_search:
            return (False, ExistanceReasonNotExisting(self))

        y = self._check_existance_find_one()
        logging.getLogger(__name__).debug('check_existance: Found one: %s', y)
        return (True, ExistanceReasonFoundOne(self, y))

    def _check_existance_primitives_more_equal_not_normals_approx(self):
        logging.getLogger(__name__).debug('_check_existance_primitives_more_equal_not_normals_approx')
        n = self.n
        q = self.q
        e = self.e
        qn = q**n

        # Check approximation for euler_phi:
        # euler_phi(n) >= n/( e^gamma * loglog n + 3/loglog n )
        euler_phi_approx = (qn-1)/(e**euler_gamma * log(log(qn-1)) + 3/log(log(qn-1)))
        if euler_phi_approx > self.ff_extension.upper_border_not_normals():
            return True

    def _factor(self):
        logging.getLogger(__name__).debug('_factor')
        # Factor q**n-1 = prod_(d|e*n) Phi_d(p)
        n = self.n
        p = self.p
        e = self.e
        return factor_with_euler_phi(p, e*n)

    def _check_existance_primitives_more_equal_not_normals(self):
        logging.getLogger(__name__).debug('_check_existance_primitives_more_equal_not_normals')
        n = self.n
        q = self.q
        e = self.e
        qn = q**n

        phi = prod((euler_phi(p)**r for p,r in self.factorization))
        if phi > self.ff_extension.upper_border_not_normals():
            return True

    def _check_existance_find_one(self):
        return self.ff_extension.pcn_element(self.factorization)

    def __str__(self):
        return "<PCNExistenceChecker q=%d^%d n=%d>" % (self.p, self.e, self.n)


if __name__ == '__main__':
    import sys
    from ff_pcn.datastore import datastore
    logging.basicConfig(level=logging.INFO)
    PCNExistenceChecker.check_range(int(sys.argv[1]), int(sys.argv[2]))
