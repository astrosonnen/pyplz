import numpy as np
from scipy.integrate import quad


Mpc = 3.08568025e24
c = 2.99792458e10
H0 = 70.
omegaL = 0.7
omegaM = 0.3

def comovd(z):
    I = quad(lambda z: 1./(omegaL + omegaM*(1+z)**3)**0.5, 0., z)
    return c/(H0*10.**5)*I[0]


