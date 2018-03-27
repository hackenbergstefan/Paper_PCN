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


YAFU_WORK_FOLDER = './yafu_job'
"""Path to working folder."""

YAFU_EXECUTABLE = './yafu-setup-package/prefix/bin/yafu'
"""Path to yafu setup package"""

TIMEOUT = 10*60


def factor_with_yafu(num,
                     timeout=TIMEOUT,
                     yafu_executable=YAFU_EXECUTABLE,
                     factor_append_to=None,
                     abort_append_to=None):

    tmpdir = tempfile.TemporaryDirectory()

    logging.critical('Start: %d', num)

    cmd = [
        os.path.abspath(yafu_executable),
        'factor(%d)' % num,
        '-of',
        'out.txt',
        '-silent'
    ]
    logging.getLogger(__name__).debug('Popen %s', ' '.join(cmd))
    proc = subprocess.Popen(
        cmd,
        cwd=tmpdir.name,
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
        out = open(tmpdir.name+'/out.txt').read()
        logging.critical('Finished: %s', out)
        if factor_append_to:
            with open(factor_append_to, 'a') as fp:
                fp.write(out+'\n')
    finally:
        tmpdir.cleanup()


def factor_with_yafu_mult(num):
    factor_with_yafu(num, factor_append_to='out', abort_append_to='abort')


def factor_batch_with_yafu(batch):
    nums = [int(n) for n in open(batch).readlines()]
    pool = multiprocessing.Pool(multiprocessing.cpu_count())

    pool.map(factor_with_yafu_mult, nums)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    factor_batch_with_yafu(sys.argv[1])
