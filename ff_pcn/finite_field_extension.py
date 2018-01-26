#!/usr/bin/env sage

"""
Module holding classes and methods modeling extensions of finite fields.
"""

__author__ = "Stefan Hackenberg"


import logging
from sage.all import moebius, divisors, Integer
from ff_pcn.finite_field_theory import essential_divisors, euler_polynomial


class FiniteFieldExtension(object):

    def __init__(self, p, r, n):
        """
        :param p: Prime p. Characteristic of Field.
        :param r: Power of base field F = GF(q), q = p^r
        :param n: Grade of extension E = GF(q^n)
        """
        self.p = Integer(p)
        self.r = Integer(r)
        self.n = Integer(n)
        self.q = p**r


    def upper_border_not_normals(self):
        """
        Returns a (very good) upper border for the number of elements of E, which are not normal.
        """
        logging.getLogger(__name__).debug('upper_border_not_normals_1')
        p = self.p
        r = self.r
        q = self.q
        n = self.n

        essential_divs = essential_divisors(p, r, n)
        border = sum(
            sum(
                moebius(n//(d*a))*q**(d*a) - euler_polynomial(q, d, n)
                for a in divisors(n//d)
            )
            for d in essential_divs
        )
        logging.getLogger(__name__).debug('upper_border_not_normals %d', border)
        return border
