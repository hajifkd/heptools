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
    eq_a(kB.in_(GeV), 1.)

def test_subunits():
    eq_(str(MeV), '1.00e+00 MeV')
    eq_(str(keV), '1.00e+00 keV')
    eq_(str(keV**2), '1.00e+00 keV^2')
    eq_(str(b), '1.00e+00 b')
    eq_a((3.89379 * 10**11 * GeV**2 * fb).in_GeV(), 1.)

def test_conv():
    eq_(cm.in_m(), 1e-3 * m)
    eq_(MeV.in_GeV(), 1e-3 * GeV)
