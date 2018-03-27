#!/usr/bin/env python2

try:
    import ff_pcn
except ImportError:
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(__file__, '../../')))

import sys
import os
import shutil
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
                     yafu_work_folder=YAFU_WORK_FOLDER,
                     factor_append_to=None,
                     abort_append_to=None):

    if os.path.exists(yafu_work_folder):
        shutil.rmtree(yafu_work_folder)
    os.mkdir(yafu_work_folder)

    logging.critical('Start: %d', num)

    cmd = [
        os.path.abspath(yafu_executable),
        'factor(%d)' % num,
        '-of',
        'out.txt',
        '-threads',
        '%d' % multiprocessing.cpu_count(),
        '-silent'
    ]
    logging.getLogger(__name__).debug('Popen %s', ' '.join(cmd))
    proc = subprocess.Popen(
        cmd,
        cwd=yafu_work_folder,
    )

    try:
        proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        logging.critical('Abort')
        if abort_append_to:
            with open(abort_append_to, 'a') as fp:
                fp.write('%d\n' % num)
    else:
        out = open(yafu_work_folder+'/out.txt').read()
        logging.critical('Finished: %s', out)
        if factor_append_to:
            with open(factor_append_to, 'a') as fp:
                fp.write(out+'\n')


def factor_batch_with_yafu(batch):
    nums = [int(n) for n in open(batch).readlines()]
    for num in nums:
        factor_with_yafu(num, factor_append_to=batch+'.out', abort_append_to=batch+'.abort')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    factor_batch_with_yafu(sys.argv[1])
