#!/usr/bin/env python2

"""
"""

__author__ = "Stefan Hackenberg"


import logging
import os
import sys
import re
import csv
from sage.all import factor, Integer, prod


FACTOR_DATABASE = os.path.abspath(os.path.join(__file__, '../factors.csv'))


def facprod(fac):
    return prod((p**m for p, m in fac))


class Factorer(object):

    def __init__(self):
        self.database = dict()
        self.load()
        self.queue = []

    def add(self, num, factorization=None):
        self.database[num] = factorization
        if factorization is not None:
            self.save((num, factorization))

    def get(self, num):
        if num < 1e10:
            return list(factor(num))
        if num in self.database:
            return self.database[num]
        logging.getLogger(__name__).critical('Factorer.get: Factorization needed : %d', num)
        return None

    def save(self, data=None):
        if data is not None:
            with open(FACTOR_DATABASE, 'a') as fp:
                writer = csv.writer(fp)
                writer.writerow(data)
            return
        with open(FACTOR_DATABASE, 'w') as fp:
            writer = csv.writer(fp)
            writer.writerows(sorted(self.database.items()))

    def load(self):
        with open(FACTOR_DATABASE, 'r') as fp:
            reader = csv.reader(fp)
            for num, fac in reader:
                try:
                    num = Integer(num)
                    fac = eval(fac)
                    assert num == facprod(fac)
                    self.database[num] = fac
                except SyntaxError:
                    logging.critical('Unable to eval: %s %s', num, fac)
                    raise
        logging.getLogger(__name__).debug('Factorer.load: Loaded %s', self.database)

    def read(self, yafu_out_fil):
        """
        Read yafu output file. Line format: (NUMBER)/FAC1/FAC2/...
        """
        for line in open(yafu_out_fil).readlines():
            num = int(re.search(r'^\((\d+)\)', line).group(1))
            facs = [(Integer(r), Integer(m or 1)) for r, m in re.findall(r'/(\d+)\^?(\d+)?', line)]

            num2 = facprod(facs)
            assert num == num2, '%d != %d = prod(%s)' % (num, num2, facs)

            facs = cleanup_factorization([(f, 1) for f in facs])
            self.add(num, facs)


factorer = Factorer()


def cleanup_factorization(factorization):
    """
    Cleanup factorization.
    """
    facs = {}
    for p, mul in factorization:
        if p in facs:
            facs[p] += mul
        else:
            facs[p] = mul
    return list(facs.items())


if __name__ == '__main__':
    factorer.read(sys.argv[1])
