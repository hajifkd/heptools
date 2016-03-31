import math


def eta_to_theta(eta):
    return 2 * math.atan(math.exp(-eta))


def theta_to_eta(theta):
    return -math.log(math.tan(theta / 2.))

