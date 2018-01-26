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
from sage.all import Integer, is_prime, euler_gamma, log, PolynomialRing, divisors, factor
from ff_pcn.basic_number_theory import is_regular
from ff_pcn.finite_field_extension import FiniteFieldExtension


class ExistanceReason(object):

    def __init__(self, checker):
        self.checker = checker

    def __repr__(self):
        return 'PCN exists for (%d, %d, %d)' % (self.checker.p, self.checker.r, self.checker.n)


class ExistanceReasonRegular(ExistanceReason):

    def __repr__(self):
        return '%s because (p,r,n) is regular' % ExistanceReason.__repr__(self)


class ExistanceReasonQBiggerN(ExistanceReason):

    def __repr__(self):
        return '%s because q > n' % ExistanceReason.__repr__(self)


class ExistanceReasonPrimitivesMoreEqualNotNormalsApprox(ExistanceReason):

    def __repr__(self):
        return '%s |P| >= |H| by approximation' % ExistanceReason.__repr__(self)


class ExistanceReasonPrimitivesMoreEqualNotNormals(ExistanceReason):

    def __repr__(self):
        return '%s |P| >= |H|' % ExistanceReason.__repr__(self)

class ExistanceReasonFoundOne(ExistanceReason):

    def __init__(self, checker, f):
        super(ExistanceReasonFoundOne, self).__init__(checker)
        self.f = f

    def __repr__(self):
        return '%s found %s' % (ExistanceReason.__repr__(self), self.f)

class ExistanceReasonNotExisting(ExistanceReason):

    def __repr__(self):
        return 'NO PCN exists for (%d, %d, %d)' % (self.checker.p, self.checker.r, self.checker.n)

class ExistanceReasonNeedFactorization(ExistanceReason):

    def __init__(self, checker, needed_factors):
        super(ExistanceReasonNeedFactorization, self).__init__(checker)
        self.needed_factors = needed_factors
        if not isinstance(needed_factors, list):
            needed_factors = [needed_factors]
        for fac in needed_factors:
            factorlib.add(fac, None)

    def __repr__(self):
        return 'For (%d, %d, %d) factorization is needed: %s' % (self.checker.p, self.checker.r, self.checker.n, self.needed_factors)


class PCNExistenceChecker(object):
    """
    Class for checking existence of PCN element in
    E/F where F = GF(q^n) and F = GF(q) with q = p^r.
    """

    @staticmethod
    def check_until_n(n):
        """
        Checks exitance of PCNs for all FiniteField extensions of degree n.
        """
        results = dict()
        for p in xrange(n//2, n):
            for e in xrange(1,n):
                if not is_prime(p):
                    continue
                if p**e > n:
                    continue
                checker = PCNExistenceChecker(p, e, n)
                res = results[(p,e,n)] = checker.check_existance()
                logging.getLogger(__name__).info('check_until_n of (%d, %d, %d) => %s', p, e, n, res)
        return results

    def __init__(self, p, r, n):
        self.p = Integer(p)
        self.r = Integer(r)
        self.n = Integer(n)
        self.q = p**r
        self.ff_extension = FiniteFieldExtension(self.p, self.r, self.n)

    def check_existance(self):
        logging.getLogger(__name__).debug('check_existance of %s', self)
        p = self.p
        q = self.q
        n = self.n
        r = self.r

        # Existance if (p,r,n) is regular
        if is_regular(p, r, 1, n, 1):
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

        result = self._check_existance_primitives_more_equal_not_normals()
        if result is True:
            logging.getLogger(__name__).debug('check_existance: |P| >= |H|')
            return (True, ExistanceReasonPrimitivesMoreEqualNotNormals(self))

        result = self._check_existance_find_one()
        if result is not None:
            logging.getLogger(__name__).debug('check_existance: Explicitly found %s', result)
            return True, ExistanceReasonFoundOne(self, result)

        return False, ExistanceReasonNotExisting(self)

    def _check_existance_primitives_more_equal_not_normals_approx(self):
        logging.getLogger(__name__).debug('_check_existance_primitives_more_equal_not_normals_approx')
        n = self.n
        q = self.q
        r = self.r
        qn = q*n

        # Check approximation for euler_phi:
        # euler_phi(n) >= n/( e^gamma * loglog n + 3/loglog n )
        euler_phi_approx = (qn-1)/(r**euler_gamma * log(log(qn-1)) + 3/log(log(qn-1)))
        if euler_phi_approx > self.ff_extension.upper_border_not_normals():
            return True

    def _check_existance_primitives_more_equal_not_normals(self):
        """
        Check |P| >= |H| by calculating |P|
        """
        r = self.r
        n = self.n
        p = self.p
        Zx = PolynomialRing(sage.all.ZZ,'x')
        factors_left = False
        factors = []
        for d in divisors(r*n):
            phi = Zx.cyclotomic_polynomial(d)(p)
            fac = factorlib.get(phi)
            if fac is not None:
                factors += fac
            else:
                logging.getLogger(__name__).debug('_check_existance_primitives_more_equal_not_normals: need factorization of: %d', phi)
                factorlib.add(phi, None)
                factors_left = True
        if factors_left:
            return False
        prims = 1
        for (f, mul) in factors:
            prims *= f**(mul-1)*(f-1)
            if prims > self.ff_extension.upper_border_not_normals_1():
                return True
        return False

    def _check_existance_find_one(self):
        """
        Prove existance of PCN element by finding one.
        """
        if factorlib.get(self.q**self.n - 1) is not None:
            return self.ff_extension.any_pcn()


    def __str__(self):
        return "<PCNExistenceChecker q=%d^%d n=%d>" % (self.p, self.r, self.n)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    PCNExistenceChecker.check_until_n(int(sys.argv[1]))
    # Factorer.process_factorlib()
