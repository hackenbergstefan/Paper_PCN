#!/usr/bin/env sage

"""
Package applying theoretical results about extensions of finite fields to
give computational results on existence of primitive complete normal basis generators.
"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(__file__, '../../')))


class ExistanceReason(object):

    def __init__(self, checker):
        self.checker = checker

    def __repr__(self):
        return 'PCN exists for (%d, %d, %d)' % (self.checker.p, self.checker.e, self.checker.n)


class ExistanceReasonRegular(ExistanceReason):

    def __repr__(self):
        return '%s because (p,e,n) is regular' % ExistanceReason.__repr__(self)


class ExistanceReasonQBiggerN(ExistanceReason):

    def __repr__(self):
        return '%s because q > n' % ExistanceReason.__repr__(self)


class ExistanceReasonPrimitivesMoreEqualNotNormalsApprox(ExistanceReason):

    def __repr__(self):
        return '%s |P| >= |H| by approximation' % ExistanceReason.__repr__(self)


class ExistanceReasonPrimitivesMoreEqualNotNormals(ExistanceReason):

    def __repr__(self):
        return '%s |P| >= |H|' % ExistanceReason.__repr__(self)

class ExistanceReasonFoundOne(ExistanceReason):

    def __init__(self, checker, f):
        super(ExistanceReasonFoundOne, self).__init__(checker)
        self.f = f

    def __repr__(self):
        return '%s found %s' % (ExistanceReason.__repr__(self), self.f)

class ExistanceReasonNotExisting(ExistanceReason):

    def __repr__(self):
        return 'NO PCN exists for (%d, %d, %d)' % (self.checker.p, self.checker.e, self.checker.n)

class ExistanceReasonNeedFactorization(ExistanceReason):

    def __init__(self, checker):
        super(ExistanceReasonNeedFactorization, self).__init__(checker)

    def __repr__(self):
        return 'For (%d, %d, %d) factorization is needed' % (self.checker.p, self.checker.e, self.checker.n)
