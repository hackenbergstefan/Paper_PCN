#!/usr/bin/env python2

"""
"""

__author__ = "Stefan Hackenberg"


import logging
import os
import csv
from sage.all import factor, Integer


FACTOR_DATABASE = os.path.abspath(os.path.join(__file__, '../factors.csv'))


class Factorer(object):

    def __init__(self):
        self.database = dict()
        self.load()

    def add(self, num, factorization=None):
        self.database[num] = factorization
        self.save()

    def get(self, num):
        if num < 1e10:
            return list(factor(num))
        if num in self.database:
            return self.database[num]
        logging.getLogger(__name__).critical('Factorer.get: Factor ecm: %d', num)
        facs = list(factor(num, method='ecm'))
        self.add(num, facs)
        logging.getLogger(__name__).critical('Factorer.get: Added %d = %s', num, facs)
        return facs

    def save(self):
        with open(FACTOR_DATABASE, 'w') as fp:
            writer = csv.writer(fp)
            writer.writerows(self.database.items())

    def load(self):
        with open(FACTOR_DATABASE, 'r') as fp:
            reader = csv.reader(fp)
            for num, fac in reader:
                try:
                    self.database[Integer(num)] = eval(fac)
                except SyntaxError:
                    logging.getLogger(__name__).critical('Unable to eval: %s %s', num, fac)
                    raise
        logging.getLogger(__name__).debug('Factorer.load: Loaded %s', self.database)


factorer = Factorer()


def cleanup_factorization(factorization):
    """
    Cleanup factorization.
    """
    facs = {}
    for p, mul in factorization:
        if p in facs:
            facs[p] = mul
        else:
            facs[p] += mul
    return list(facs.items())
