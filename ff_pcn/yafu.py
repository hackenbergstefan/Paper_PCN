#!/usr/bin/env python2

try:
    import ff_pcn
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(__file__, '../../')))

import sys
import os
import tempfile
import subprocess
import logging
import multiprocessing
import argparse


YAFU_WORK_FOLDER = './yafu_job'
"""Path to working folder."""

YAFU_EXECUTABLE = './yafu-setup-package/prefix/bin/yafu'
"""Path to yafu setup package"""

TIMEOUT = 60

YAFU_ARGS = []


def factor_with_yafu(num,
                     timeout=None,
                     yafu_executable=YAFU_EXECUTABLE,
                     factor_append_to=None,
                     abort_append_to=None):

    timeout = timeout or TIMEOUT

    with tempfile.TemporaryDirectory() as tmpdir:
        logging.critical('Start: %d with timeout %d', num, timeout)

        cmd = [
            os.path.abspath(yafu_executable),
            'factor(%d)' % num,
            '-of',
            'out.txt',
        ] + YAFU_ARGS

        logging.getLogger(__name__).debug('Popen %s', ' '.join(cmd))
        proc = subprocess.Popen(
            cmd,
            cwd=tmpdir,
        )

        try:
            proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.communicate()
            logging.critical('Abort %d', num)
            if abort_append_to:
                with open(abort_append_to, 'a') as fp:
                    fp.write('%d\n' % num)
        else:
            if os.path.exists(tmpdir+'/out.txt'):
                out = open(tmpdir+'/out.txt').read()
            else:
                out = '({0})/{0}'.format(num)
            logging.critical('Finished: %s', out)
            if factor_append_to:
                with open(factor_append_to, 'a') as fp:
                    fp.write(out+'\n')


def factor_with_yafu_mult(num):
    factor_with_yafu(num, factor_append_to='out', abort_append_to='abort')


def factor_batch_with_yafu(batch, threads):
    nums = [int(n) for n in open(batch).readlines()]
    pool = multiprocessing.Pool(threads)

    pool.map(factor_with_yafu_mult, nums)


def main():
    global TIMEOUT
    global YAFU_ARGS
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument(
        '--timeout',
        type=int,
        default=TIMEOUT,
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=multiprocessing.cpu_count(),
    )
    parser.add_argument(
        '--yafu-args',
        default='--silent',
    )

    args = parser.parse_args()
    TIMEOUT = args.timeout
    YAFU_ARGS = args.yafu_args.split()
    factor_batch_with_yafu(args.file, threads=args.threads)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
