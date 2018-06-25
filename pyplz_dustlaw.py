#this program encodes the dust law from Cardelli et al. 1989
#takes a wavelength (in Angstrom) and a redshift and spits out the ratio between (rest frame) A(lambda) and A(V)

xfunc = lambda l: 10000./l
yfunc = lambda l: xfunc(l)-1.82
dfunc = lambda p,l: p[0] + p[1]*yfunc(l) + p[2]*yfunc(l)**2 + p[3]*yfunc(l)**3 + p[4]*yfunc(l)**4 + p[5]*yfunc(l)**5 + p[6]*yfunc(l)**6 + p[7]*yfunc(l)**7
pol7 = lambda p,y: p[0] + p[1]*y + p[2]*y**2 + p[3]*y**3 + p[4]*y**4 + p[5]*y**5 + p[6]*y**6 + p[7]*y**7

#values of the coefficient of the dust law function from Cardelli et al. 1989
pa = [1, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999]
pb = [0, 1.41338, 2.2830, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002]


def Fa(x):
    if x <= 8. and x >= 5.9:
        return -0.04473*(x-5.9)**2 + 0.1207*(x-5.9)**3
    elif x < 5.9:
        return 0.
def Fb(x):
    if x <= 8. and x >= 5.9:
        return 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3
    elif x < 5.9:
        return 0.

def a(x):
    if x >=0.3 and x <= 1.1: #infrared
        return 0.574*x**1.61
    elif x > 1.1 and x <= 3.3: #Optical/near UV
        y = x - 1.82
        return pol7(pa,y)
    elif x > 3.3 and x <= 8.: #UV
        return 1.752 - 0.316*x - 0.104/((x-4.67)**2+0.341) + Fa(x)

def b(x):
    if x >=0.3 and x <= 1.1: #infrared
        return -0.527*x**1.61
    elif x > 1.1 and x <= 3.3: #Optical/near UV
        y = x - 1.82
        return pol7(pb,y)
    elif x > 3.3 and x <= 8.: #UV
        return -3.090 + 1.825*x + 1.206/((x-4.62)**2 + 0.263) + Fb(x)

#Alambda gives the ratio between the extinction at wavelength lambda A(lambda) and the extinction in the V band (rest frame)
#z is the redshift of the dust screen (default:z=0,galactic)
def Alambda(l,z=0.0, Rv=3.1):
    lrest = l/(1.+z)
    x = xfunc(lrest)
    return a(x) + b(x)/Rv

