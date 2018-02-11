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


RESULT_FOLDER = os.path.abspath(os.path.join(__file__, '../../result/'))


class Database(object):

    re_filename = re.compile(r'ex_(?P<n>\d*)\.txt')

    def __init__(self, result_folder=RESULT_FOLDER):
        self.result_folder = result_folder

    def add(self, p, e, n, result):
        with open(os.path.join(self.result_folder, 'ex_%d.txt' % n), 'a') as fp:
            fp.write('%s\n' % result)

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

    def find_missing_pcns(self, fil):
        from ff_pcn.pcn_existence_checker import PCNExistenceChecker
        from sage.all import Integer
        with open(fil, 'r') as fp:
            content = fp.readlines()

        for i, line in enumerate(content):
            if 'False' in line:
                break
        if i == len(content):
            return

        logging.getLogger(__name__).info('check line: %s' % line)
        match = re.search(r'^\((\d+), (\d+), (\d+)\)', line)
        p, e, n = tuple(map(Integer, match.groups()))
        checker = PCNExistenceChecker(p, e, p**e, n)
        result = checker.check_existance(no_explicit_search=False)
        logging.getLogger(__name__).info(result)
        with open(fil, 'w') as fp:
            fp.writelines(content[:i] + ['%s\n' % result[1]] + content[i+1:])


database = Database()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    # database.missing_files(int(sys.argv[1]))
    # database.check_and_cleanup()
    database.find_missing_pcns(sys.argv[1])
