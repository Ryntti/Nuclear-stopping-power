
import numpy as np
import constants as c
from math import sqrt
import conversions as conv
from scipy.optimize import fsolve
from numpy.polynomial.legendre import leggauss 

# screening length function
def scr_length(z1: int, z2: int):
    return 0.46848*10**(-10)/(z1**0.23 + z2**0.23)

# the screening function
def scr_function(x: float):
    """
    The screening function. Only called in this module when setting up the
    potential function. Important to note is that the units are in SI.
    """

    sum = 0
    for i in range(len(c.alpha)):
        sum += c.alpha[i]*np.exp(-c.beta[i]*x)
    return sum

# The potential V
def V(r: float, z1: int, z2: int):
    """
    Screened Coulomb potential. The first argument is the radial distance
    from the source, the second and third are the atomic mass numbers of the
    colliding particles. 
    Returns in SI units i.e. J/C

    """

    scr_arg = r/scr_length(z1, z2)
    return scr_function(scr_arg)*(z1*z2*c.q**2)/(4*np.pi*c.e0*r)

# the function g(r) squared. It's only used for determining r_min.
def g_squared(r: float, b: float, E_com: float, z1: int, z2: int):
    """
    The function g(r) squared. g(r) is only directly called when solving for
    rmin in the equation g(r) = 0. Since g(r) is a square root of a rational 
    function, g(r) and the rational function have the same roots. In fact, the 
    rational function is just g squared.

    Arguments assumed in SI units, notably that E_com is in Joules.
    """

    return 1 - (b/r)**2 - V(r, z1, z2)/E_com 


def F(u: float, b: float, E_com: float, r_min: float, z1: int, z2: int):
    """  
    The integrand of the scattering angle formula. u has been substituted in 
    place of r to help the integration.

    Arguments assumed in SI units, notably that E_com is in Joules.
    """

    return 1/sqrt(b**2*(2 - u**2) + r_min**2/(u**2*E_com)*(V(r_min, z1, z2) - V(r_min/(1-u**2), z1, z2)))


