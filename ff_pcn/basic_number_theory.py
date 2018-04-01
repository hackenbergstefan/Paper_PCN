#!/usr/bin/env python2

"""
Module coding basic number theoretical applications.
"""

__author__ = "Stefan Hackenberg"


from sage.all import (
    GF,
    Hom,
    Integer,
    cyclotomic_polynomial,
    divisors,
    factor,
    gcd,
    is_prime,
    moebius,
    prime_divisors,
    prod,
    uniq,
)
from ff_pcn import MissingFactorsException


def regular(p, e, n):
    """
    Returns True if (q, n) with q = p^e is regular.
    """
    p = Integer(p)
    e = Integer(e)
    n = Integer(n)
    return is_regular(p, e, 1, n, 1)


def is_regular(p, e, k, t, pi):
    """
    Tests if (p,e,k,t,pi) is regular
    """
    return gcd(ordn(squarefree(k*p_free_part(t, p)), p**e), k*t*pi) == 1


def squarefree(n):
    """
    Returns squarefree part of n. Also called nu(n).
    """
    return prod(map(lambda x: x[0], factor(Integer(n))))


def ordn(m, q):
    """
    Computes ordn_m(q) = min{ k: q ** k = 1 mod m }
    """
    if m == 1 or q == 1:
        return 1

    q_ = q % m
    for i in xrange(1, m+1):
        if q_ == 1:
            return i
        q_ = (q_ * q) % m


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


def factor_with_euler_phi(p, m, use_factorer=True):
    """
    Returns factorization of p**m-1 with p prime by using
    p**m-1 = prod_(d|m) Phi_d(p).
    """
    from ff_pcn.factorer import factorer
    pm = p**m-1
    factors = []
    missing_factors = []
    for d in divisors(m):
        phi = cyclotomic_polynomial(d)(p)
        assert phi.divides(pm)
        if phi == 1:
            continue
        if use_factorer:
            facs = factorer.get((d, p))
        else:
            facs = list(factor(phi))
        factors += facs or []
        if facs is None:
            missing_factors += [(d, p, phi)]
    if len(missing_factors) == 0:
        ret = {}
        for k, l in factors:
            if k not in ret:
                ret[k] = l
            else:
                ret[k] += l
        assert prod(k**l for k, l in ret.items()) == p**m - 1
        return sorted(ret.items())
    else:
        factorer.queue += missing_factors
    raise MissingFactorsException(missing_factors)


def euler_phi(factorization):
    """
    Returns euler_phi(n) by giving a factorization of n.
    """
    return prod(
        p**(k-1) * (p-1)
        for p, k
        in factorization
    )


def cyclotomic_equivalents(n, b):
    """
    Returns equivalent pairs (n',b') for (n,b) by using identity

    Phi_mp^l(b) = Phi_mp(b^(p^(l-1))).
    (n,b) is always included.
    """
    ret = [(n, b)]
    qs = [q**(m-1) for q, m in list(factor(n)) if m > 1]
    for q in qs:
        ret += cyclotomic_equivalents(n//q, b**q)
    return sorted(uniq(ret), key=lambda nb: (nb[1], nb[0]))
