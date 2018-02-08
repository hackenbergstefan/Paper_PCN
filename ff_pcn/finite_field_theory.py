#!/usr/bin/env sage

"""
Module coding results mainly found by Hachenberger.
"""

__author__ = "Stefan Hackenberg"


import logging
import itertools
from sage.all import Integer, uniq, factor, DiGraph, is_prime, divisors, prime_divisors, prod, euler_phi
from ff_pcn.basic_number_theory import largest_divisor, ordn, squarefree, p_free_part, is_regular
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


@store('euler_polynomial')
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

    return q**(n-d*tau) * \
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
    if is_regular(p, e, 1, n, 1):
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
