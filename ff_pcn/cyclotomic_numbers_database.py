#!/usr/bin/env python


"""
Module downloading factorizations from
http://www.asahi-net.or.jp/~KC2H-MSM/cn/
"""

import re
import requests
from sage.all import (
    Integer,
    euler_phi,
    is_prime,
    cyclotomic_polynomial,
    prod,
)


DATABSE_URL = "http://www.asahi-net.or.jp/~KC2H-MSM/cn/"


re_database_line = re.compile(r'\((?P<n>\d+) (?P<b>\d+) \((?P<fac>.+)\) \(P \d+\)\)')


def get_factorization(n, b):
    """
    Returns factorization of Phi_n(b) from database. If not existing None is returned.
    """
    n = Integer(n)
    b = Integer(b)
    phin = euler_phi(n)
    req = requests.get(DATABSE_URL + 'p%dn%03d.txt' % (phin, n))
    if req.status_code != 200:
        req = requests.get(DATABSE_URL + 'old/p%dn%03d.txt' % (phin, n))
        if req.status_code != 200:
            return None
    fac = None
    for line in req.content.splitlines():
        match = re_database_line.match(line)
        if match and int(match.group('b')) == b:
            fac = [Integer(p) for p in match.group('fac').split()]
            break
    if fac is None:
        return
    phi = cyclotomic_polynomial(n)(b)
    fac += [phi//prod(fac)]
    fac = [p for p in fac if p != 1]
    assert phi == prod(fac)
    assert all(map(is_prime, fac))
    return [(p, 1) for p in fac]
