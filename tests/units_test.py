from heptools.units import *
from nose.tools import ok_, eq_
from . import eq_a

def test_add():
    eq_(GeV(1) + GeV(1.5), GeV(2.5))

def test_div():
    eq_(GeV(5) / GeV(1), 5.)

def test_mult():
    eq_(3 * GeV(1), GeV(3))

def test_consts():
    eq_a((0.197 * GeV * 10**-15 * m).in_(GeV), 1.)
    eq_a(c.in_(GeV), 1.)
    eq_a(kb.in_(GeV), 1.)
