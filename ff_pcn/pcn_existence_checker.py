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
from sage.all import Integer, euler_gamma, log, PolynomialRing, divisors, factor, primes, ZZ, prod, euler_phi, uniq
from ff_pcn import ExistanceReasonRegular, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox, ExistanceReasonPrimitivesMoreEqualNotNormals, ExistanceReasonNeedFactorization, ExistanceReasonFoundOne, ExistanceReasonNotExisting, MissingFactorsException, ExistanceReasonProposition53
from ff_pcn.basic_number_theory import is_regular, factor_with_euler_phi, p_free_part
from ff_pcn.finite_field_extension import FiniteFieldExtension
from ff_pcn.database import database
from ff_pcn.factorer import factorer
from ff_pcn.finite_field_theory import pens_to_check


def check_p_n(pn):
    p, n = pn
    for e in xrange(1,n):
        q = p**e
        if q >= p_free_part(n, p):
            continue
        if q >= n:
            break
        if is_regular(p, e, 1, n, 1):
            continue
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
        for n in xrange(l, m):
            check_n(n)
        # pool = multiprocessing.Pool(multiprocessing.cpu_count())
        # pool.map(check_n, xrange(l, m))

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
                    if q >= p_free_part(n, p):
                        continue
                    if q >= n:
                        break
                    if is_regular(p, e, 1, n, 1):
                        continue
                    checker = PCNExistenceChecker(p, e, q, n)
                    res = results[(p,e,n)] = checker.check_existance()
                    logging.getLogger(__name__).info('check_to (%d, %d, %d) => %s', p, e, n, res)
        logging.getLogger(__name__).info('Missing: %d', len([r for r in results.values() if isinstance(r[1], ExistanceReasonNotExisting)]))
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
        result = self._proposition53_i()
        if result is True:
            logging.getLogger(__name__).debug('check_existance: Proposition 5.3 ')
            return (True, ExistanceReasonProposition53(self))

        result = self._proposition53_ii()
        if result is True:
            logging.getLogger(__name__).debug('check_existance: Proposition 5.3 ')
            return (True, ExistanceReasonProposition53(self))

        result = self._check_existance_primitives_more_equal_not_normals_approx()
        if result is True:
            logging.getLogger(__name__).debug('check_existance: |P| >= |H| approx')
            return (True, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox(self))

        if no_explicit_search:
            return (False, ExistanceReasonNotExisting(self))

        ret = self.factorization = self._factor()
        if isinstance(ret, tuple) and ret[0] is False:
            self.missing_factors = ret[1]
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
        p = self.p
        e = self.e
        qn = p**(e*n)

        # Check approximation for euler_phi:
        # euler_phi(n) >= n/( e^gamma * loglog n + 3/loglog n )
        euler_phi_approx = (qn-1)/(e**euler_gamma * log(log(qn-1)) + 3/log(log(qn-1)))
        if euler_phi_approx > u_qn(p, e, n):
            return True

    def _proposition53_i(self):
        logging.getLogger(__name__).debug('_proposition53')
        n = self.n
        p = self.p
        e = self.e
        q = self.q
        ls = q**n - u_qn(p, e, n)
        rs = 4514.7 * q**(5*n/8) * 2**sum(
            omega_d(d, p, e, n)
            for d in
            essential_divisors(p, e, n)
        )
        return ls >= rs

    def _proposition53_ii(self):
        logging.getLogger(__name__).debug('_proposition53')
        n = self.n
        p = self.p
        e = self.e
        q = self.q
        ls = q**n - u_qn(p, e, n)
        rs = 4514.7 * q**(5*n/8) * prod(
            theta_d(d, p, e, n) * 2 ** omega_d(d, p, e, n)
            for d in
            essential_divisors(p, e, n)
        )
        return ls >= rs

    def _factor(self):
        logging.getLogger(__name__).debug('_factor')
        # Factor q**n-1 = prod_(d|e*n) Phi_d(p)
        n = self.n
        p = self.p
        e = self.e
        try:
            return factor_with_euler_phi(p, e*n)
        except MissingFactorsException as e:
            return False, e.missing_factors

    def _check_existance_primitives_more_equal_not_normals(self):
        logging.getLogger(__name__).debug('_check_existance_primitives_more_equal_not_normals')
        n = self.n
        q = self.q
        e = self.e

        phi = prod((l**(k-1)*(l-1) for l, k in self.factorization))
        if phi > u_qn(self.p, e, n):
            return True

    def _check_existance_find_one(self):
        return self.ff_extension.pcn_element(self.factorization)

    def __str__(self):
        return "<PCNExistenceChecker q=%d^%d n=%d>" % (self.p, self.e, self.n)


class CriterionChecker(object):

    def __init__(self, pens):
        for p, e, n in pens:
            self.check_criterions(p, e, n)

    def check_criterions(self, p, e, n):
        ff = FiniteFieldExtension(p, e, n)
        crits = [
            ff.pcn_criterion_1(),
            ff.pcn_criterion_2(),
            ff.pcn_criterion_3(),
        ]
        if not any(crits):
            try:
                crits += [
                    ff.pcn_criterion_4(),
                ]
            except MissingFactorsException:
                crits += [None]
        logging.info('check_criterions %s: %s', (p, e, n), crits)
        if not any(crits):
            logging.critical('check_criterions %s: None True', (p, e, n))

    def any_criterion(self, p, e, n):
        ff = FiniteFieldExtension(p, e, n)
        ret = any((
            ff.pcn_criterion_1(),
            ff.pcn_criterion_2(),
            ff.pcn_criterion_3(),
            ff.pcn_criterion_4(),
        ))
        if not ret:
            logging.info('any_criterion %s', (p, e, n))
        return ret


if __name__ == '__main__':
    import sys
    # from ff_pcn.datastore import datastore
    logging.basicConfig(level=logging.INFO)
    # PCNExistenceChecker.check_range(int(sys.argv[1]), int(sys.argv[2]))
    # PCNExistenceChecker.check_to(int(sys.argv[1]))
    for n in xrange(int(sys.argv[1]), int(sys.argv[2])):
        pens = pens_to_check(n)
        CriterionChecker(pens)

    queue = ['factor(%d)' % n for n in sorted(uniq(factorer.queue))]
    logging.critical('factorizations needed: \n%s', '\n'.join(queue))
