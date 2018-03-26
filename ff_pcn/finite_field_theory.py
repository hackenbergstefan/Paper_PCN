#!/usr/bin/env python2

"""
Module coding results mainly found by Hachenberger.
"""

__author__ = "Stefan Hackenberg"


import logging
import itertools
from sage.all import (
    DiGraph,
    Integer,
    divisors,
    euler_gamma,
    euler_phi,
    e as euler_const,
    factor,
    is_prime,
    log,
    moebius,
    prime_divisors,
    prod,
    primes,
    uniq,
)
from ff_pcn.basic_number_theory import largest_divisor, ordn, squarefree, p_free_part, regular
from ff_pcn.datastore import store
from ff_pcn.factorer import factorer


def decompose(p, e, n):
    """
    Application of the Decomposition Theorem (Section 19) for x^n-1 over F_p^e.
    """
    pi = largest_divisor(p, n)
    return decompose_cyclic_module(p, e, 1, n/pi, pi)


def decompose_cyclic_module(p, e, k, t, pi):
    """
    Internal application of the Decomposition Theorem for Phi_k(x^(t*pi)) over F_p^e.
    """
    logging.getLogger(__name__).debug('decompose_cyclic_module (%d, %d, (%d,%d,%d))', p, e, k, t, pi)
    assert not p.divides(k*t), 'p must not divide kt'

    for r, l in reversed(factor(t)):
        if not (r**l).divides(ordn(squarefree(k*t), p**e)):
            R = largest_divisor(r, t)
            return decompose_cyclic_module(p, e, k, t/r, pi) + decompose_cyclic_module(p, e, k*R, t/R, pi)
    return [(k, t, pi)]


def module_characters(decomp):
    """
    Returns the module characters of a given decomposition:
    The module character of U_F,Phi_k(x^t) is k*t / nu(k)
    """
    return uniq(map(lambda l: l[0]*l[1]*l[2] / squarefree(l[0]), decomp))


# @store('euler_polynomial')
def euler_polynomial(q, d, n):
    """
    Returns phi_(q^d)(x^(n/d)-1) where phi_q is the polynomial analogon for the euler totient function.

    :param q: Cardinality of a finite field (i.e. must be a prime power).
    :param d: Degree of extension of GF(q).
    :param n: Degree of total extension.
    """
    p = prime_divisors(q)[0]
    tau = p_free_part(n//d, p)

    qd = q**d
    divs_tau = divisors(tau)
    ordne = {e: ordn(e, qd) for e in divs_tau}

    return q**(n//d-tau) * \
        prod((
              (qd**ordne[e] - 1)**(euler_phi(e)//ordne[e])
              for e in divisors(tau)
          ))


def universal_essential_set(n):
    """
    Returns the universal essential set of divisors of n.
    """
    n = Integer(n)
    prims = prime_divisors(n)
    prims_good = filter(lambda r: not any([r.divides(s-1) for s in prims]), prims)
    prims_good = dict((r,get_multiplicity(r,n)) for r in prims_good)
    divsN = divisors(n)[:-1]
    adjfunc = (lambda i,j:
               i.divides(j) and (Integer(j/i) in prims_good) and
               ((get_multiplicity(Integer(j/i),i) == prims_good[Integer(j/i)]-1
                 and get_multiplicity(Integer(j/i),j) == prims_good[Integer(j/i)])
                or
                (get_multiplicity(Integer(j/i),i) == prims_good[Integer(j/i)]-2
                 and get_multiplicity(Integer(j/i),j) == prims_good[Integer(j/i)]-1)))
    g = DiGraph([divsN, adjfunc])
    verts_indegzero = filter(lambda v: g.in_degree(vertices=[v]) == [0], g.vertices())
    return verts_indegzero


def essential_divisors(p, e, n):
    """
    Returns a list of essential divisors of (p, e, n).
    """
    q = p**e
    if regular(p, e, n):
        return []
    divsN = divisors(n)[:-1]
    adjfunc = (lambda i,j:
               i.divides(j) and is_prime(Integer(j/i)) and
               not Integer(j/i).divides(ordn(p_free_part(n/j, p), q**i)))
    g = DiGraph([divsN, adjfunc])
    verts_indegzero = filter(lambda v: g.in_degree(vertices=[v]) == [0], g.vertices())
    divsModChar = list(uniq(itertools.chain(*map(divisors, module_characters(decompose(p, e, n))))))
    essential_divs = filter(lambda d: d in divsModChar, verts_indegzero)
    logging.getLogger(__name__).debug('essential_divisors (%d, %d, %d) => %s', p, e, n, essential_divs)
    return essential_divs


def primitive_element(E, facs):
    """
    Returns a primitive element of E, for given factorization of |E|-1
    """
    qn = E.order() - 1
    cofacs = [qn//p for p, mul in facs]
    for y in E:
        if y in [0, 1]:
            continue
        if all((y**co) != 1 for co in cofacs):
            return y


def u_qn(p, e, n):
    """
    Returns U_(p**e,n). Proposition 4.4.
    """
    q = p**e

    essential_divs = essential_divisors(p, e, n)
    border = sum(
        sum(
            moebius(n//(d*a))*q**(d*a)
            for a in divisors(n//d)
        ) - euler_polynomial(q, d, n)
        for d in essential_divs
    )
    assert border >= 0
    return border


def omega_d(d, p, e, n):
    """
    Returns Omega_d := sum_(t|(n/d)') phi(t)/ord_t(q^d).
    """
    return sum([
        euler_phi(t)//ordn(t, p**(e*d))
        for t in
        divisors(p_free_part(n//d, p))
    ])


def theta_d(d, p, e, n):
    """
    Returns Theta_d := Phi_(q^d)(x^(n/d)' - 1) / q^(d*(n/d)').
    """
    q = p**e
    n_ = p_free_part(n//d, p)*d
    ret = euler_polynomial(q, d, n_) / q**(d * p_free_part(n//d, p))
    assert ret < 1
    return ret


def lower_euler_phi(n):
    """
    Returns a lower bound for euler_phi(n).
    """
    return n/(euler_const**euler_gamma * log(log(n)) + 3/log(log(n)))


def pens_to_check(n):
    """
    Returns a list of tripples (p, e, n) with:
      - p**e < n'
      - (p, e, n) is not regular
    """
    n = Integer(n)
    tocheck = []
    for p in primes(n):
        for e in xrange(1, n):
            if regular(p, e, n):
                continue
            if p**e >= p_free_part(n, p):
                continue
            tocheck += [(p, e, n)]
    return tocheck
