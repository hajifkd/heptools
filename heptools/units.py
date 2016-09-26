#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import functools
import sys

from sympy import Rational

# this might be useful
from math import pi

BASE = 'GeV'


class UnitNotMatchError(Exception):
    pass


def sameunits(f):
    @functools.wraps(f)
    def wrapper(u1, u2):
        if not u1.has_same_units(u2):
            raise UnitNotMatchError()
        return f(u1, u2)
    return wrapper


class Unit(object):

    def __init__(self, coeff, **kwargs):
        self.coeff = coeff
        self.units = {}
        for k in kwargs:
            if kwargs[k]:
                self.units[k] = kwargs[k]

    def __new__(cls, coeff, **kwargs):
        if any(kwargs.values()):
            return super(Unit, cls).__new__(cls)
        else:
            return coeff

    def __str__(self):
        res = '%3.2e' % self.coeff
        for k in self.units:
            if self.units[k] == 1:
                res += ' %s' % k
            else:
                res += ' %s^%s' % (k, str(self.units[k]))
        return res

    def __repr__(self):
        res = '%3.2e' % self.coeff
        for k in self.units:
            if self.units[k] == 1:
                res += ' * %s' % k
            else:
                res += ' * %s**(%s)' % (k, str(self.units[k]))
        return res

    def __call__(self, coeff):
        return coeff * self

    def has_same_units(self, t):
        unit1 = set(self.units)
        unit2 = set(t.units)

        if unit1 != unit2:
            return False

        for u in unit1:
            if self.units[u] != t.units[u]:
                return False

        return True

    def to_natural_unit(self):
        p = 0
        coeff = self.coeff
        for k in self.units:
            if k == 'GeV':
                p += self.units[k]
            else:
                cc, cp = Unit.unity[k]
                coeff *= cc**self.units[k]
                p += cp * self.units[k]

        return coeff, p

    def convert(self, t):
        sc, sp = self.to_natural_unit()
        tc, tp = t.to_natural_unit()

        return (sc / tc**(1.0 * sp / tp), Rational(sp, tp))

    def in_(self, t):
        c, p = self.convert(t)
        return c * t**p

    def __getattr__(self, name):
        if name.startswith('in_'):
            u = getattr(sys.modules[__name__], name[3:])
            return functools.partial(self.in_, u)
        else:
            raise AttributeError

    def inverse(self):
        return 1.0 / self

    @sameunits
    def __add__(self, t):
        return Unit(self.coeff + t.coeff, **self.units)

    @sameunits
    def __sub__(self, t):
        return Unit(self.coeff - t.coeff, **self.units)

    def __pow__(self, t):
        return Unit(self.coeff**t,
                    **{k: self.units[k] * t for k in self.units})

    def __mul__(self, t):
        if isinstance(t, Unit):
            u = {}
            for k in self.units:
                u[k] = self.units[k]
            for k in t.units:
                if k in u:
                    u[k] += t.units[k]
                else:
                    u[k] = t.units[k]
            return Unit(self.coeff * t.coeff, **u)
        else:
            return Unit(self.coeff * t, **self.units)

    def __rmul__(self, t):
        return self * t

    def __div__(self, t):
        if isinstance(t, Unit):
            u = {}
            for k in self.units:
                u[k] = self.units[k]
            for k in t.units:
                if k in u:
                    u[k] -= t.units[k]
                else:
                    u[k] = -t.units[k]

            return Unit(self.coeff / t.coeff, **u)
        else:
            return Unit(self.coeff / t, **self.units)

    __truediv__ = __div__

    def __rdiv__(self, t):
        return Unit(t / self.coeff, **{k: -self.units[k] for k in self.units})

    __rtruediv__ = __rdiv__

    @sameunits
    def __cmp__(self, t):
        return self.coeff - t.coeff

    def __hash__(self):
        return self.coeff.__hash__() ^ self.units.__hash__()

# units
GeV = Unit(1, GeV=1)
MeV = Unit(1e-3, GeV=1)
keV = Unit(1e-6, GeV=1)

m = Unit(1, m=1)
km = Unit(1e3, m=1)
cm = Unit(1e-2, m=1)

s = Unit(1, s=1)
kg = Unit(1, kg=1)
kelvin = Unit(1, kelvin=1)

# physical constants
c = 299792458 * m / s
hbar = 1.0545718e-34 * m**2 * kg / s
hbarc = 1.97326979e-16 * m * GeV
kb = 8.61733034e-14 * GeV / kelvin
G = 6.67408e-11 * m**3 * kg**-1 * s**-2


def __generate_unity():
    unity = {}
    # m
    unity['m'] = (1.0 / hbarc.coeff, -hbarc.units[BASE])
    # s
    hbarGeV = hbarc / c  # s GeV
    unity['s'] = (1.0 / hbarGeV.coeff, -hbarGeV.units[BASE])
    # kg
    he = c * hbarc / hbar  # GeV/kg
    unity['kg'] = (he.coeff, he.units[BASE])
    # kelvin
    unity['kelvin'] = (kb.coeff, kb.units[BASE])

    return unity

Unit.unity = __generate_unity()

# additional constants
Mp = (8 * pi * G).in_(GeV).inverse()**Rational(1, 2)
