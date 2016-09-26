from nose.tools import ok_

def eq_a(a, b, u=1):
    ok_(a - b < 0.01 * u and a - b > -0.01 * u)
