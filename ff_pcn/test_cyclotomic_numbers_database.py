#!/usr/bin/env python2

"""
Test for basic_number_theory.
"""

from unittest import TestCase
from ff_pcn.cyclotomic_numbers_database import get_factorization


class CyclotomicNumbersDatabaseTestCase(TestCase):

    def test_get_factorization(self):
        for n, b in [(61, 50), (220, 59)]:
            self.assertGreater(len(get_factorization(n, b)), 0)
