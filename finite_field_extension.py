#!/usr/bin/env sage

"""
Module holding class representing a FiniteFieldExtension.
"""

__author__ = "Stefan Hackenberg"


import logging
import itertools
import sage.all
from sage.all import gcd, Integer, factor, divisors, prime_divisors, uniq, moebius, GF, PolynomialRing, Hom, is_prime, euler_phi

from utils import squarefree, ordn, largest_divisor, is_regular

def get_multiplicity(p,n):
    """
    returns multiplicity of p in n
    """
    a = 0
    while p.divides(n):
        a += 1
        n = sage.Integer(n/p)
    return a

def get_module_characters(decomp):
    """
    returns the module characters of a decomposition
    the module character of U_F,Phi_k(x^t) is k*t / nu(k)
    """
    return uniq(map(lambda l: l[0]*l[1]*l[2] / squarefree(l[0]),decomp))

def decompose(p,e, n):
    """
    Application of the Decomposition Theorem (Section 19)
    for x^n-1 over F_p^e
    """
    pi = largest_divisor(p,n)
    return decompose_cycl_module(p,e, 1, n/pi, pi)

def decompose_cycl_module(p,e, k,t,pi):
    """
    internal application of the Decomposition Theorem
    for Phi_k(x^(t*pi)) over F_p^e
    """
    logging.getLogger(__name__).debug('decompose_cycl_module (%d, %d, (%d,%d,%d))', p, e, k, t, pi)
    if p.divides(k*t):
        print "ERROR p | kt"
    #test all prime divisors, start with largest one
    flag = False
    for r,l in reversed(factor(t)):
        if not (r**l).divides(ordn(squarefree(k*t),p**e)):
            R = largest_divisor(r,t)
            return decompose_cycl_module(p,e, k, t/r, pi) \
                    + decompose_cycl_module(p,e, k*R, t/R, pi)
    return [(k,t,pi)]

def get_universal_essential_set(n):
    n = sage.Integer(n)
    prims = prime_divisors(n)
    prims_good = filter(lambda r: not any([r.divides(s-1) for s in prims]), prims)
    prims_good = dict((r,get_multiplicity(r,n)) for r in prims_good)
    #print "prims_good", prims_good
    divsN = divisors(n)[:-1]
    adjfunc = (lambda i,j:
            i.divides(j) and (Integer(j/i) in prims_good) and\
            ((get_multiplicity(Integer(j/i),i) == prims_good[Integer(j/i)]-1\
             and get_multiplicity(Integer(j/i),j) == prims_good[Integer(j/i)])\
            or
            (get_multiplicity(Integer(j/i),i) == prims_good[Integer(j/i)]-2\
             and get_multiplicity(Integer(j/i),j) == prims_good[Integer(j/i)]-1)\
            ))
    g = sage.all.DiGraph([divsN, adjfunc])
    verts_indegzero = filter(lambda v: \
            g.in_degree(vertices=[v]) == [0], g.vertices())
    #g.show()
    return verts_indegzero


def essential_divisors(p,e,n):
    logging.getLogger(__name__).debug('essential_divisors (%d, %d, %d)', p, e, n)
    p = Integer(p)
    e = Integer(e)
    n = Integer(n)
    q = Integer(p**e)
    if is_regular(p,e,1,n,1): return "isRegular"
    divsN = divisors(n)[:-1]
    adjfunc = (lambda i,j:\
        i.divides(j) and is_prime(Integer(j/i)) and \
                not Integer(j/i).divides(ordn(p_free_part(n/j,p),q**i)))
    g = sage.all.DiGraph([divsN, adjfunc])
    #g.show()
    verts_indegzero = filter(lambda v: \
            g.in_degree(vertices=[v]) == [0], g.vertices())
    #adjfunc2 = (lambda j,i:\
        #sage.Integer(i).divides(j) and is_prime(sage.Integer(j/i)) and \
        #ordn(p_free_part(n/i,p),q**i) == sage.Integer(j/i)*ordn(p_free_part(n/j,p),q**i) )
    #g2 = sage.DiGraph([verts_indegzero,adjfunc2])
    #g2.show()
    #verts_indegzero2 = filter(lambda v: \
            #g2.in_degree(vertices=[v]) == [0], g2.vertices())
    #if verts_indegzero2 != verts_indegzero:
        #print "case found on: ",(p,e,n)
        #return g,g2
    divsModChar = list(uniq(itertools.chain(*map(divisors,\
        get_module_characters(decompose(p,e,n))))))
    essential_divs = filter(lambda d: d in divsModChar, verts_indegzero)
    logging.getLogger(__name__).debug('essential_divisors (%d, %d, %d) => %s', p, e, n, essential_divs)
    return essential_divs


class NotDecomposableException(Exception):
    pass


def decompose_cyclic_module(p,e, k,t,pi):
    """
    internal application of the Decomposition Theorem
    for Phi_k(x^(t*pi)) over F_p^e
    """
    if p.divides(k*t):
        raise NotDecomposableException('p divides k*t')

    # Test all prime divisors, start with largest one
    flag = False
    for r, l in reversed(factor(t)):
        if not (r**l).divides(ordn(squarefree(k*t),p**e)):
            R = largest_divisor(r,t)
            return decompose_cyclic_module(p,e, k, t/r, pi) + decompose_cyclic_module(p,e, k*R, t/R, pi)
    return [(k,t,pi)]


def p_free_part(t,p):
    """
    Returns p-free part of t
    """
    while p.divides(t):
        t /= p
    return t


def not_completely_basic_divisors(p,e,n):
    """
    Returns the NOT completely basic divisors of an extension n over GF(p^e)
    """
    n = Integer(n)
    q = Integer(p**e)
    divs = []
    divsN = divisors(n)
    while len(divsN) > 0:
        d = divsN.pop(0)
        is_compl_basic = True
        for r in prime_divisors(n/d):
            if r.divides(ordn(p_free_part(n/d/r,p),q**d)):
                is_compl_basic = False
                break
        divs += [d]
        if is_compl_basic:
            divsN = filter(lambda k: not d.divides(k), divsN)
    return divs


def module_characters(decomp):
    """
    Returns the module characters of a decomposition.
    The module character of U_F,Phi_k(x^t) is k*t / nu(k)
    """
    return uniq(map(lambda l: l[0]*l[1]*l[2] / squarefree(l[0]), decomp))


def euler_polynomial(q, d, n):
    """
    returns phi_(q^d)(x^(n/d)-1)
    where phi_q is the polynomial analogon for the euler totient function

    q: cardinality of a finite field (i.e. must be a prime power)
    d: degree of extension of GF(q)
    n: degree of total extension
    """
    if not hasattr(euler_polynomial, 'cache'):
        setattr(euler_polynomial, 'cache', dict())
    cached = euler_polynomial.cache.get((q, d, n), None)
    if cached is not None:
        return cached

    logging.getLogger(__name__).debug('euler_polynomial (%d, %d, %d)', q, d, n)
    p = prime_divisors(q)[0]
    tau = Integer(n/d)
    while p.divides(tau): tau = Integer(tau/p)

    ret = q**(n-d*tau) * \
            sage.all.prod( [ (q**(d*ordn(e,q**d)) - 1)**(euler_phi(e)/ordn(e,q**d))\
            for e in divisors(tau) ] )
    logging.getLogger(__name__).debug('euler_polynomial (%d, %d, %d) => %s', q, d, n, ret)

    euler_polynomial.cache[(q, d, n)] = ret
    return ret


def completely_normal(x, E, q, divs, fieldsAll, facsAll, prodsAll):
    """
    Tests if x in E is CN over GF(q).
    """
    if x == E.zero():
        return False
    #test isNormal for each divisor
    pows = dict()
    for d in divs:
        h = Hom(fieldsAll[d],E)[0];
        for idx,(f,mult) in enumerate(facsAll[d]):
            g = prodsAll[d][idx];
            ret = E.zero();
            iold = 0
            xiold = x
            for i,gi in enumerate(list(g)):
                if pows.has_key(i*d):
                    xi = pows[i*d];
                    iold = i*d
                    xiold = xi
                else:
                    xi = xiold**(q**(d*i-iold));
                    pows[i*d] = xi;
                    xiold = xi
                    iold = i*d
                ret += h(gi)*xi
            if ret == 0: return False;
    return True


class FiniteFieldExtension(object):

    def __init__(self, p, r, n):
        """
        :param p: Prime p. Characteristic of Field.
        :param r: Power of base field F = GF(q), q = p^r
        :param n: Grade of extension E = GF(q^n)
        """
        self.p = Integer(p)
        self.r = Integer(r)
        self.n = Integer(n)
        self.q = p**r

    def proper_subfield_divisors(self):
        divs = getattr(self, '_cache_proper_subfield_divisors', None)
        if divs is not None:
            return divs
        p = self.p
        e = self.r
        n = self.n
        basList = not_completely_basic_divisors(p,e,n)
        divs = filter(
            lambda d: d in basList,
            list(uniq(itertools.chain(
                *map(divisors, module_characters(self.decompose()))))
            )
        )
        setattr(self, '_cache_proper_subfield_divisors', divs)
        return divs

    def decompose(self):
        """
        Application of Decomposition Theorem (Hachenberger ...)
        for x^n-1 over F = GF(p^r)
        """
        p = self.p
        r = self.r
        n = self.n
        pi = largest_divisor(p, n)
        return decompose_cyclic_module(p, r, 1, n/pi, pi)


    def upper_border_not_normals(self):
        """
        Returns an upper border for the number of elements of E, which are not normal.
        """
        logging.getLogger(__name__).debug('upper_border_not_normals_1')
        p = self.p
        r = self.r
        q = self.q
        n = self.n
        essential_divs = essential_divisors(p, r, n)

        logging.getLogger(__name__).debug('%s', [[n//(d*a) for a in divisors(n//d)] for d in essential_divs])
        ret = sum(
            sum(
                moebius(n//(d*a))*q**(d*a) - euler_polynomial(q, d, n)
                for a in divisors(n//d)
            )
            for d in essential_divs
        )
        logging.getLogger(__name__).debug('ret %d', ret)
        return ret

    def any_pcn(self, primitive_element=None):
        """
        Searches for a PCN element in E|F.
        """
        q = self.q
        n = self.n

        # Setup required fields
        F = GF(q)
        E = F.extension(n,'a')
        P = E.prime_subfield()

        Px = PolynomialRing(P,'x')
        Fx = PolynomialRing(F,'x')
        Ex = PolynomialRing(E,'x')
        h = Hom(F,E)[0]
        prim_order = E.order()-1
        primitives = []

        # Setup factors of x^n-1
        divs = self.proper_subfield_divisors()
        facs_all = dict();
        prods_all = dict();
        fields_all = dict();
        for d in divs:
            G = F.extension(Integer(d), 'c');
            Gx = PolynomialRing(G,'x');
            fields_all[d] = G;
            facs_all[d] = list((Gx.gen()**(n/d)-1).factor());
            prods_all[d] = dict();
            for idx,(f,mult) in enumerate(facs_all[d]):
                prods_all[d][idx] = (Gx.gen()**(n/d)-1).quo_rem(f)[0]

        # get one primitive element
        if primitive_element == None:
            x = E.primitive_element()
        else:
            x = primitive_element
            E = primitive_element.parent()

        lasti = 0
        y = x
        for i in itertools.count(1):
            if gcd(i, E.order()-1) != 1: continue
            y = y*x**(i-lasti)
            lasti = i
            if completely_normal(y, E, q, divs, fields_all, facs_all, prods_all):
                mipo = y.minpoly()
                for f,i in Fx(mipo).factor():
                    if f.map_coefficients(h)(y) == E.zero():
                        return f
