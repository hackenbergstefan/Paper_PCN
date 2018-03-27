#! /usr/bin/env python

import sys
import os
import subprocess
import logging
import multiprocessing


YAFU_WORK_FOLDER = './yafu_job'

TIMEOUT = 10*60


def factor_with_yafu(num, timeout=TIMEOUT):
    with open('job.bat', 'w') as fp:
        fp.write('factor(%d)\n' % num)

    proc = subprocess.Popen(
        ['yafu', '-batchfile', 'job.bat', '-of', 'out.txt', '-threads', '%d' % multiprocessing.cpu_count()],
        cwd=YAFU_WORK_FOLDER,
        stdout=subprocess.PIPE,
    )

    try:
        proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
    else:
        logging.info('Finished: %s', open(YAFU_WORK_FOLDER+'/out.txt').readlines())


if __name__ == '__main__':
    factor_with_yafu(int(sys.argv[1]))
