#!/usr/bin/env python2

"""
Test for basic_number_theory.
"""

from unittest import TestCase
from ff_pcn.basic_number_theory import (
    squarefree,
    ordn,
    p_free_part,
    regular,
    factor_with_euler_phi,
    euler_phi as euler_phi_from_factor,
    cyclotomic_equivalents,
)
from sage.all import (
    lcm,
    euler_phi,
    Integer,
    primes,
    factor,
    cyclotomic_polynomial
)


class BasicNumberTheoryTestCase(TestCase):

    def test_squarefree(self):
        for x, y in [(1, 1), (4, 2), (24, 6)]:
            self.assertEqual(squarefree(x), y)

    def test_ordn(self):
        self.assertEqual(ordn(1, 1), 1)
        self.assertEqual(ordn(2, 3), 1)
        self.assertEqual(ordn(8*17, 7**2), lcm(ordn(8, 7**2), ordn(17, 7**2)))
        self.assertTrue(Integer(ordn(8*17, 7**2)).divides(euler_phi(8*17)))

    def test_p_free_part(self):
        p = Integer(17)
        self.assertEqual(p_free_part(p, p), 1)
        self.assertFalse(p.divides(p_free_part(p**3*17, p)))
        self.assertEqual(p_free_part(Integer(23), p), 23)

    def test_regular(self):
        self.assertFalse(regular(2, 3, 12))
        for n in [15, 33, 35, 45, 51, 65, 69, 75, 77, 85, 87, 91, 95, 99, 115, 119, 123, 133, 135, 141, 143, 145, 153, 159, 161, 175, 177, 185, 187]:
            for p in primes(100):
                for e in xrange(5):
                    self.assertTrue(regular(p, e, n), str((p,e,n)))

    def test_factor_with_euler_phi(self):
        for p in primes(50):
            for e in xrange(1, 5):
                self.assertEqual(factor_with_euler_phi(p, e, use_factorer=False), list(factor(p**e-1)))

    def test_euler_phi(self):
        for n in xrange(1, 10000):
            n = Integer(n)
            self.assertEqual(euler_phi_from_factor(list(factor(n))), euler_phi(n))

    def test_cyclotomic_equivalents(self):
        n = Integer(2**3 * 3**2)
        b = Integer(5)
        equivs = cyclotomic_equivalents(n, b)
        self.assertEqual(equivs, [(n, b), (n//3, b**3), (n//4, b**4), (n//12, b**12)])
        self.assertEqual(len(set([cyclotomic_polynomial(n)(b) for n, b in equivs])), 1)
