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
        if not any(crits) and None not in crits:
            crits += [
                ff.pcn_criterion_5(),
            ]
        if not any(crits) and None not in crits:
            crits += [
                ff.pcn_criterion_6(),
            ]
        logging.info('check_criterions %s: %s', (p, e, n), crits)
        if not any(crits):
            logging.critical('check_criterions %s: None True', (p, e, n))

if __name__ == '__main__':
    import sys
    # from ff_pcn.datastore import datastore
    logging.basicConfig(level=logging.INFO)
    # PCNExistenceChecker.check_range(int(sys.argv[1]), int(sys.argv[2]))
    # PCNExistenceChecker.check_to(int(sys.argv[1]))
    for n in xrange(int(sys.argv[1]), int(sys.argv[2])):
        pens = pens_to_check(n)
        CriterionChecker(pens)

    queue = ['%d %d %d %d' % (euler_phi(d), d, p, phi) for d, p, phi in sorted(uniq(factorer.queue), key=lambda (d, p, phi): euler_phi(d))]
    if len(queue):
        logging.critical('factorizations needed: \n%s', '\n'.join(queue))
