#!/usr/bin/env python2

"""
"""

__author__ = "Stefan Hackenberg"


try:
    import ff_pcn
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(__file__, '../../')))
import logging
import os
import sys
import re
import csv
from sage.all import factor, Integer, prod, cyclotomic_polynomial
from ff_pcn.cyclotomic_numbers_database import get_factorization as get_factorization_from_online_database
from ff_pcn.basic_number_theory import cyclotomic_equivalents


FACTOR_DATABASE = os.path.abspath(os.path.join(__file__, '../cyclotomic_numbers.csv'))


def facprod(fac):
    return prod((p**m for p, m in fac))


class Factorer(object):

    def __init__(self):
        self.database = dict()
        self.load()
        self.queue = []

    def add(self, nb, num, factorization):
        self.database[nb] = factorization
        self.save()

    def get(self, nb):
        num = cyclotomic_polynomial(nb[0])(nb[1])
        if num < 1e10:
            return list(factor(num))

        equivalents = cyclotomic_equivalents(*nb)

        # Lookup local database
        for mb in equivalents:
            if mb in self.database:
                if nb != mb:
                    self.database[nb] = self.database[mb]
                    del self.database[mb]
                return self.database[nb]

        # Lookup online database
        for mb in equivalents:
            fac = get_factorization_from_online_database(*mb)
            logging.getLogger(__name__).critical('Factorer.get: online lookup: %s %s', mb, fac)
            if fac is not None:
                self.add(nb, num, fac)
                return fac

        logging.getLogger(__name__).critical('Factorer.get: Factorization needed : %s %d', nb, num)
        return None

    def save(self):
        with open(FACTOR_DATABASE, 'w') as fp:
            writer = csv.writer(fp)
            writer.writerows(sorted(self.database.items()))

    def load(self):
        with open(FACTOR_DATABASE, 'r') as fp:
            reader = csv.reader(fp)
            for nb, fac in reader:
                try:
                    nb = eval(nb)
                    fac = eval(fac)
                    num = cyclotomic_polynomial(nb[0])(nb[1])
                    assert num == facprod(fac)
                    self.database[nb] = fac
                except SyntaxError:
                    logging.critical('Unable to eval: %s %s', nb, fac)
                    raise
        logging.getLogger(__name__).debug('Factorer.load: Loaded %s', self.database)

    def read(self, yafu_out_fil):
        """
        Read yafu output file. Line format: (NUMBER)/FAC1/FAC2/...
        """
        for line in open(yafu_out_fil).readlines():
            match = re.search(r'^\((\d+), (\d+)\) \((\d+)', line)
            if not match:
                continue
            n = int(match.group(1))
            b = int(match.group(2))
            num = int(match.group(3))
            facs = [(Integer(r), Integer(m or 1)) for r, m in re.findall(r'/(\d+)\^?(\d+)?', line)]

            num2 = facprod(facs)
            assert num == num2, '%d != %d = prod(%s)' % (num, num2, facs)

            facs = cleanup_factorization(facs)
            self.add((n, b), num, facs)
        self.save()


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
