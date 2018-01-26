#!/usr/bin/env sage

"""
Module containing helper functionality.
"""

__author__ = "Stefan Hackenberg"


import sage.all
from sage.all import gcd, Integer


def is_regular(p,e, k,t,pi):
    return gcd( ordn( squarefree(k*p_free_part(t,p)), p**e ),  k*t*pi) == 1


def p_free_part(t,p):
    while p.divides(t):
        t /= p
    return t


def ordn(m,q):
    """
    Computes ordn m(q) = min{ k: q ** k = 1 mod m }
    """
    if m == 1: return 1
    for i in range(1,m+1):
        if (q ** i)%m == 1: return i;


def squarefree(n):
    return sage.all.prod(map(lambda x: x[0], sage.all.factor(Integer(n))))


def notnormals(p, r, n):
    divs = get_proper_subfield_divisors(p,r,n)
    notnorms = sum([n/k*(q**(n-k)-1) for k in divs])


def largest_divisor(p, n):
    """
    Returns largest power of p dividing n.
    """
    l = 0
    while (p**l).divides(n):
        l = l+1
    return p**(l-1);
