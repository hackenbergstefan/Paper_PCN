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
    sqrt,
)
from ff_pcn.basic_number_theory import (
    regular,
    factor_with_euler_phi,
    euler_phi,
)
from ff_pcn.finite_field_theory import (
    essential_divisors,
    lower_euler_phi,
    omega_d,
    primitive_element,
    is_primitive,
    theta_d,
    u_qn,
    completely_normal,
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

    def pcn_element(self):
        """
        Returns a PCN-Element.

        :param facs: Factorization of q**n-1.
        """
        self._setup_pcn_search()

        q = self.q
        n = self.n
        order = q**n - 1

        x = primitive_element(self.E, self.factorization)

        y = self.E(1)
        for i in itertools.count(1):
            y *= x
            if gcd(i, order) != 1:
                continue
            if self.completely_normal(y):
                return y

    def pcn_polynom(self):
        """
        Returns lexicographic smallest polynom in F[x] of degree n with pcn root.
        """
        def polynom_candidates(fx, deg):
            return (
                fx.gen()**deg + a * fx.gen()**(deg-1) + f
                for a, f in
                itertools.product(fx.base_ring(), fx.polynomials(max_degree=deg - 2))
                if a != 0
            )
        fx = PolynomialRing(GF(self.p, 'a'), 'x')
        logging.getLogger(__name__).debug('pcn_polynom')
        for f in polynom_candidates(fx, self.e*self.n):
            if not f.is_irreducible():
                continue
            logging.getLogger(__name__).debug('pcn_polynom: test f = %s', f)
            y = GF(self.q**self.n, name='a', modulus=f).gen()
            if is_primitive(y, self.factorization) and completely_normal(self.p, self.e, self.n, y):
                return f

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
        if hasattr(self, 'factorization'):
            return self.factorization
        self.factorization = factor_with_euler_phi(self.p, self.e*self.n, use_factorer=use_factorer)
        return self.factorization

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
        rs = sqrt(self.q**self.n).n() * (2.0**omega - 1) * prod(
            self.theta_d(d) * 2.0 ** self.omega_d(d)
            for d in
            self.essential_divisors()
        )
        logging.getLogger(__name__).debug('pcn_criterion_4: %E >= %E', ls, rs)
        assert rs >= 0
        return ls >= rs

    def pcn_criterion_5(self):
        """
        Returns True, if Criterion 5 applies.

        NOTE: Factorization of q^n-1 required.

        Criterion 5:
        phi(q^n-1) > U_qn
        """
        ls = euler_phi(self.factorization)
        rs = self.u_qn()
        logging.getLogger(__name__).debug('pcn_criterion_5: %E >= %E', ls, rs)
        assert ls > 0
        assert rs > 0
        return ls > rs

    def pcn_criterion_6(self):
        """
        Returns smallest PCN polynom, if Criterion 6 applies.

        NOTE: Factorization of q^n-1 required.

        Criterion 6:
        Give explicit PCN polynom, i.e. the lexicographic smallest polynom in F[x] of degree n, with pcn root.
        """
        return self.pcn_polynom()
