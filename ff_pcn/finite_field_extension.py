#!/usr/bin/env sage

"""
Module holding classes and methods modeling extensions of finite fields.
"""

__author__ = "Stefan Hackenberg"


import itertools
import logging
from sage.all import moebius, divisors, Integer, GF, gcd, Hom, PolynomialRing, factor
from ff_pcn.finite_field_theory import essential_divisors, euler_polynomial, primitive_element


class FiniteFieldExtension(object):

    def __init__(self, p, e, q, n):
        """
        :param p: Prime p. Characteristic of Field.
        :param e: Power of base field F = GF(q)
        :param q: q = p^e
        :param n: Grade of extension E = GF(q^n)
        """
        self.p = p
        self.e = e
        self.n = n
        self.q = q

    def upper_border_not_normals(self):
        """
        Returns a (very good) upper border for the number of elements of E, which are not normal.
        """
        logging.getLogger(__name__).debug('upper_border_not_normals_1')
        p = self.p
        e = self.e
        q = self.q
        n = self.n

        essential_divs = essential_divisors(p, e, n)
        border = sum(
            sum(
                moebius(n//(d*a))*q**(d*a)
                for a in divisors(n//d)
            ) - euler_polynomial(q, d, n)
            for d in essential_divs
        )
        logging.getLogger(__name__).debug('upper_border_not_normals %d', border)
        return border

    def pcn_element(self, facs):
        """
        Returns a PCN-Element.

        :param facs: Factorization of q**n-1.
        """
        q = self.q
        n = self.n
        F = GF(q, 'a')
        E = F.extension(n, 'a')
        order = q**n - 1

        x = primitive_element(E, facs)

        y = E(1)
        for i in itertools.count(1):
            y *= x
            if gcd(i, order) != 1:
                continue
            if self.completely_normal(y):
                return y

    def completely_normal(self, y):
        """
        Returns True if y in GF(q**n) is completely normal over GF(q)
        """
        p = self.p
        e = self.e
        q = self.q
        n = self.n

        if y in [0, 1]:
            return False

        self._pow_cache = {1: y}

        divs = essential_divisors(p, e, n)
        return all(self.normal(y, d) for d in divs)

    def normal(self, y, d):
        """
        Returns True if y in E is normal over G = GF(q^d).

        Theorem: y is normal over G iff:
        f(sigma^d)(y) != 0
        for all f cofactors of x^(n/d) - 1 over G
        with sigma: x -> x^q.
        """
        logging.getLogger(__name__).debug('normal: test y: %s over extension %s over F', y, d)
        n = self.n
        G = GF(self.q).extension(d)
        qd = G.order()
        Gx = PolynomialRing(G, 'x')
        h = Hom(G, y.parent())[0]
        basepol = Gx.gen()**(n//d)-1
        cofactors = [basepol.quo_rem(f)[0] for f, mul in list(basepol.factor())]

        # Test if frobenius vanishes on cofactors
        for cofac in cofactors:
            # cofac = cofac(Gx.gen()**qd)
            logging.getLogger(__name__).debug('normal: test cofac: %s', cofac)
            if cofac.map_coefficients(h)(y**qd) == 0:
                return False
        return True

    def _cache_poly_eval(self, pol, y):
        """
        Evaluates pol(y) using cached data.
        """
        ret = y.parent().zero()
        for power, coeff in reversed(list(pol)):
            x = self._pow_cache.get(power, None)
            if x is None:
                largest_cached_pow = sorted(self._pow_cache.keys())[-1]
                x = self._pow_cache[power] = self._pow_cache[largest_cached_pow]**(power-largest_cached_pow)
            ret += coeff * x
        return ret
