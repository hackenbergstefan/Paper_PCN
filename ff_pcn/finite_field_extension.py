#!/usr/bin/env sage

"""
Module holding classes and methods modeling extensions of finite fields.
"""

__author__ = "Stefan Hackenberg"


import itertools
import logging
from sage.all import (
    GF,
    Hom,
    Integer,
    PolynomialRing,
    gcd,
    prod,
)
from ff_pcn.basic_number_theory import (
    regular,
    factor_with_euler_phi,
)
from ff_pcn.finite_field_theory import (
    essential_divisors,
    lower_euler_phi,
    omega_d,
    primitive_element,
    theta_d,
    u_qn,
)


class FiniteFieldExtension(object):

    def __init__(self, p, e, n):
        """
        :param p: Prime p. Characteristic of Field.
        :param e: Power of base field F = GF(q)
        :param q: q = p^e
        :param n: Grade of extension E = GF(q^n)
        """
        self.p = Integer(p)
        self.e = Integer(e)
        self.n = Integer(n)
        self.q = self.p**self.e
        self.qn = self.q**self.n
        self._cache_cofactors = dict()

    def pcn_element(self, facs):
        """
        Returns a PCN-Element.

        :param facs: Factorization of q**n-1.
        """
        self._setup_pcn_search()

        q = self.q
        n = self.n
        order = q**n - 1

        x = primitive_element(self.E, facs)

        y = self.E(1)
        for i in itertools.count(1):
            y *= x
            if gcd(i, order) != 1:
                continue
            if self.completely_normal(y):
                return y

    def _setup_pcn_search(self):
        """
        Setup caches for completely normal tests.
        """
        q = self.q
        p = self.p
        e = self.e
        n = self.n
        self.essential_divs = essential_divisors(p, e, n)
        self.F = GF(q, 'a')
        self.E = self.F.extension(n, 'a')
        self.cofactors = dict()
        for d in self.essential_divs:
            G = self.F.extension(d)
            Gx = PolynomialRing(G, 'x')
            h = Hom(G, self.E)[0]
            basepol = Gx.gen()**(n//d)-1
            cofactors = [basepol.quo_rem(f)[0] for f, mul in list(basepol.factor())]
            cofactors = [f.map_coefficients(h) for f in cofactors]
            self.cofactors[d] = cofactors
            del h
            del Gx
            del G

    def completely_normal(self, y):
        """
        Returns True if y in GF(q**n) is completely normal over GF(q)
        """
        if not hasattr(self, 'essential_divs'):
            self._setup_pcn_search()

        if y in [0, 1]:
            return False

        return all(self.normal(y, d) for d in self.essential_divs)

    def normal(self, y, d):
        """
        Returns True if y in E is normal over G = GF(q^d).

        Theorem: y is normal over G iff:
        f(sigma^d)(y) != 0
        for all f cofactors of x^(n/d) - 1 over G
        with sigma: x -> x^q.
        """
        logging.getLogger(__name__).debug('normal: test y: %s over extension %s over F', y, d)

        if not hasattr(self, 'cofactors'):
            self._setup_pcn_search()

        yd = y**(self.q*d)

        # Test if frobenius vanishes on cofactors
        for cofac in self.cofactors[d]:
            logging.getLogger(__name__).debug('normal: test cofac: %s', cofac)
            if cofac(yd) == 0:
                return False
        return True

    def _cache_poly_eval(self, pol, y):
        """
        Evaluates pol(y) using cached data.
        """
        ret = y.parent().zero()
        for power, coeff in enumerate(reversed(list(pol))):
            x = self._pow_cache.get(power, None)
            if x is None:
                largest_cached_pow = sorted(self._pow_cache.keys())[-1]
                x = self._pow_cache[power] = self._pow_cache[largest_cached_pow]**(power-largest_cached_pow)
            ret += coeff * x
        return ret

    def omega_d(self, d):
        """
        Returns Omega_d := sum_(t|(n/d)') phi(t)/ord_t(q^d).
        """
        assert Integer(d).divides(self.n)
        return omega_d(d, self.p, self.e, self.n)

    def theta_d(self, d):
        """
        Returns Theta_d := Phi_(q^d)(x^(n/d)' - 1) / q^(d*(n/d)').
        """
        assert Integer(d).divides(self.n)
        return theta_d(d, self.p, self.e, self.n)

    def u_qn(self):
        """
        Returns U_(p**e,n). Proposition 4.4.
        """
        return u_qn(self.p, self.e, self.n)

    def l_qn(self):
        """
        Returns L_(p**e,n). Equation 4.2.
        """
        return lower_euler_phi(self.qn - 1).n(10)

    def essential_divisors(self):
        """
        Returns a list of essential divisors of (p, e, n).
        """
        return essential_divisors(self.p, self.e, self.n)

    def regular(self):
        return regular(self.p, self.e, self.n)

    def factor(self, use_factorer=True):
        """
        Returns factorization of q^n - 1.
        """
        return factor_with_euler_phi(self.p, self.e*self.n, use_factorer=use_factorer)

    def pcn_criterion_1(self):
        """
        Returns True, if Criterion 1 applies.

        Criterion 1:
        L_qn > U_qn
        """
        ls = self.l_qn()
        rs = self.u_qn()
        logging.getLogger(__name__).debug('pcn_criterion_1: %E > %E', ls, rs)
        return ls > rs

    def pcn_criterion_2(self):
        """
        Returns True, if Criterion 2 applies.

        Criterion 2: Equation 5.3
        q^n - U_qn >= 4514.7 * q^(5n/8) 2^sum_d Omega_d
        """
        qn = self.qn
        ls = self.qn - self.u_qn()
        rs = 4514.7 * qn**(5.0/8) * 2**sum(
            self.omega_d(d)
            for d in
            self.essential_divisors()
        )
        logging.getLogger(__name__).debug('pcn_criterion_2: %E >= %E', ls, rs)
        assert rs >= 0
        return ls >= rs

    def pcn_criterion_3(self):
        """
        Returns True, if Criterion 3 applies.

        Criterion 3: Equation 5.3 with 5.1
        q^n - U_qn >= 4514.7 * q^(5n/8) prod_d( Theta_d * 2^Omega_d )
        """
        qn = self.qn
        ls = self.qn - self.u_qn()
        rs = 4514.7 * qn**(5.0/8) * prod(
            self.theta_d(d) * 2 ** self.omega_d(d)
            for d in
            self.essential_divisors()
        )
        logging.getLogger(__name__).debug('pcn_criterion_3: %E >= %E', ls, rs)
        assert rs >= 0
        if ls < rs:
            rs = 4.9 * qn**(3.0/4) * prod(
                self.theta_d(d) * 2 ** self.omega_d(d)
                for d in
                self.essential_divisors()
            )
            logging.getLogger(__name__).debug('pcn_criterion_3: %E >= %E', ls, rs)
        return ls >= rs

    def pcn_criterion_4(self):
        """
        Returns True, if Criterion 3 applies.

        NOTE: Factorization of q^n-1 required.

        Criterion 4:
        Equation 5.1 with q^n - U_qn
        """
        factorization = self.factor()
        omega = len(factorization)
        ls = self.qn - self.u_qn()
        rs = self.q**(self.n/2.0) * (2**omega - 1) * prod(
            self.theta_d(d) * 2 ** self.omega_d(d)
            for d in
            self.essential_divisors()
        )
        logging.getLogger(__name__).debug('pcn_criterion_4: %E >= %E', ls, rs)
        assert rs >= 0
        return ls >= rs
