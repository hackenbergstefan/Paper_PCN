#!/usr/bin/env python2

"""
Database storing results.
"""

__author__ = "Stefan Hackenberg"

try:
    import ff_pcn
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(__file__, '../../')))
import os
import re
import logging
from ff_pcn import ExistanceReasonRegular, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox, ExistanceReasonPrimitivesMoreEqualNotNormals, ExistanceReasonNeedFactorization, ExistanceReasonFoundOne, ExistanceReasonNotExisting


RESULT_FOLDER = os.path.abspath(os.path.join(__file__, '../../result/'))


class Database(object):

    re_filename = re.compile(r'ex_(?P<n>\d*)\.txt')

    def __init__(self, result_folder=RESULT_FOLDER):
        self.result_folder = result_folder

    def add(self, p, e, n, result):
        resstring = ''
        if isinstance(result, ExistanceReasonRegular):
            resstring = 'regular'
        elif isinstance(result, ExistanceReasonPrimitivesMoreEqualNotNormalsApprox):
            resstring = 'L > U'
        elif isinstance(result, ExistanceReasonPrimitivesMoreEqualNotNormals):
            resstring = '|P| > |H|'
        elif isinstance(result, ExistanceReasonNeedFactorization):
            resstring = 'Factorization needed'
        elif isinstance(result, ExistanceReasonFoundOne):
            resstring = 'found %s' % result
        elif isinstance(result, ExistanceReasonNotExisting):
            resstring = 'not existing'
        with open(os.path.join(self.result_folder, 'ex_%d.txt' % n), 'a') as fp:
            fp.write('(%d, %d, %d) => %s\n' % (p, e, n, resstring))

    def check_and_cleanup(self):
        for fil in sorted(
                os.listdir(self.result_folder),
                key=lambda fil: int(self.re_filename.match(fil).group('n')),
                reverse=True):
            if self.re_filename.match(fil):
                self.check_and_cleanup_file(os.path.join(self.result_folder, fil))

    def check_and_cleanup_file(self, fil):
        from ff_pcn.pcn_existence_checker import qs_to_check
        logging.info('check_and_cleanup_file %s', fil)
        n = long(self.re_filename.match(os.path.basename(fil)).group('n'))
        qs = qs_to_check(n)
        content = open(fil, 'r').read()
        lines = content.splitlines()

        if len(lines) < len(qs):
            logging.critical('too less lines %d < %d', len(lines), len(qs))

        results = []
        for p, e in qs:
            match = re.search(r'^\(%d, %d,? %d\)(.*?)$' % (p, e, n), content, re.MULTILINE)
            if not match:
                logging.critical('(%d, %d, %d) missing', p, e, n)
            logging.debug('(%d, %d, %d) => %s', p, e, n, match.group(0))
            results += [
                (p, e, '(%d, %d, %d)%s' % (p, e, n, match.group(1)))
            ]

        results = sorted(results)
        with open(fil, 'w') as fp:
            fp.write('\n'.join(r[2] for r in results))

    def missing_files(self, m):
        for n in xrange(1, m):
            if not os.path.exists(os.path.join(self.result_folder, 'ex_%d.txt' % n)):
                logging.critical('missing %d', n)


database = Database()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    database.missing_files(int(sys.argv[1]))
    database.check_and_cleanup()
