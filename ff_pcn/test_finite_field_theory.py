#!/usr/bin/env python2

"""
Test for basic_number_theory.
"""

from unittest import TestCase
from ff_pcn.finite_field_theory import (
    decompose,
    euler_polynomial,
    u_qn,
    lower_euler_phi,
)
from sage.all import (
    Integer,
    euler_phi,
    primes,
)


class FiniteFieldTheoryTestCase(TestCase):

    def test_decompose(self):
        self.assertEqual(
            decompose(Integer(3), Integer(1), Integer(20)),
            [(1, 1, 1), (2, 1, 1), (4, 1, 1), (5, 4, 1)]
        )
        self.assertEqual(
            decompose(Integer(5), Integer(1), Integer(252)),
            [(1, 1, 1), (2, 1, 1), (4, 1, 1), (3, 2, 1), (12, 1, 1), (9, 2, 1), (36, 1, 1), (7, 6, 1), (28, 3, 1), (63, 2, 1), (252, 1, 1)]
        )

    def test_euler_polynomial(self):
        self.assertEqual(euler_polynomial(Integer(2), Integer(1), Integer(6)), 24)

    def test_u_qn(self):
        self.assertGreaterEqual(
            12,
            2**6 - u_qn(Integer(2), Integer(1), Integer(6))
        )

    def test_lower_euler_phi(self):
        for p in primes(50):
            for e in xrange(5):
                for n in xrange(5):
                    self.assertGreaterEqual(
                        euler_phi(p**(e*n)-1),
                        lower_euler_phi(p**(e*n)-1)
                    )
