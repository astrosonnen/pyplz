import numpy as np
import cosmology
from scipy.interpolate import splrep, splev, splint
from scipy.integrate import quad
from sonnentools.cgsconstants import *
import pylab


Mpc = 3.08568025e24
c = 2.99792458e10
H0 = 70.
omegaL = 0.7
omegaM = 0.3

def comovd(z):
    I = quad(lambda z: 1./(omegaL + omegaM*(1+z)**3)**0.5, 0., z)
    return c/(H0*10.**5)*I[0]

# luminosity function from Arnouts et al. (2005)
z_arn = [0., 0.3, 0.5, 0.7, 1.0]
phistar_arn = [4.07, 6.15, 1.69, 1.67, 1.14]
Mstar_arn = [-18.05, -18.38, -19.49, -19.84, -20.11]
alpha_arn = [-1.21, -1.19, -1.55, -1.60, -1.63]

# luminosity function from Oesch et al. (2010)
z_oesch= [0.75, 1.25, 1.75, 2.5]
phistar_oesch = [10.**0.52, 10.**0.90, 10.**0.63, 10.**0.49]
Mstar_oesch = [-19.17, -20.08, -20.17, -20.69]
alpha_oesch = [-1.52, -1.84, -1.60, -1.73]

# luminosity function from Bowens et al. (2015)
z_bow = [3., 3.8, 4.9]
phistar_bow = [1.71, 1.97, 0.79]
Mstar_bow = [-20.97, -20.88, -21.10]
alpha_bow = [-1.73, -1.64, -1.76]

z_bins = z_arn + z_oesch[1:] + z_bow
phistar_bins = phistar_arn + phistar_oesch[1:] + phistar_bow
Mstar_bins = Mstar_arn + Mstar_oesch[1:] + Mstar_bow
alpha_bins = alpha_arn + alpha_oesch[1:] + alpha_bow

phistar_spline = splrep(z_bins, phistar_bins, k=1)
Mstar_spline = splrep(z_bins, Mstar_bins, k=1)
alpha_spline = splrep(z_bins, alpha_bins, k=1)

def phi(M, z):
    phistar_here = splev(z, phistar_spline)
    Mstar_here = splev(z, Mstar_spline)
    alpha_here = splev(z, alpha_spline)
    return 0.4*np.log(10.) * phistar_here * 10.**(0.4*(alpha_here+1)*(Mstar_here-M)) * np.exp(-10.**(0.4*(Mstar_here-M))) * 1e-3
    

