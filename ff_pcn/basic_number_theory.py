#!/usr/bin/env sage

"""
Module coding basic number theoretical applications.
"""

__author__ = "Stefan Hackenberg"


from ff_pcn.datastore import store


from sage.all import gcd, Integer, factor, divisors, prime_divisors, uniq, moebius, GF, PolynomialRing, Hom, is_prime, euler_phi, prod


def is_regular(p, e, k, t, pi):
    """
    Tests if (p,e,k,t,pi) is regular
    """
    return gcd(ordn(squarefree(k*p_free_part(t,p)), p**e), k*t*pi) == 1


def squarefree(n):
    """
    Returns squarefree part of n. Also called nu(n).
    """
    return prod(map(lambda x: x[0], factor(Integer(n))))


@store('ordn')
def ordn(m, q):
    """
    Computes ordn m(q) = min{ k: q ** k = 1 mod m }
    """
    if m == 1:
        return 1
    for i in xrange(1, m+1):
        if (q**i) % m == 1:
            return i


def p_free_part(t, p):
    """
    Computes the p-free part of t.
    """
    while p.divides(t):
        t /= p
    return t


def multiplicity(p, n):
    """
    Returns multiplicity of p in n
    """
    a = 0
    while p.divides(n):
        a += 1
        n = Integer(n/p)
    return a


def largest_divisor(p, n):
    """
    Returns largest power of p dividing n.
    """
    return p**multiplicity(p, n)
