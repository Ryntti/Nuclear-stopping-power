import numpy as np
import conversions as conv
import helper_functions as hf
from scipy.optimize import fsolve

# the stopping power integral
def sp(abscissas: list, weights: list, E_lab: float, z1: int, z2: int, m1: float, m2: float, gamma: float):
    """  
    Calculates the stopping power integral using the Gauss-Legendre quadrature method. 
    Nodes and weights are inserted to the function as the first two arguments as lists.
    Initial projectile energy in lab coordinates, atomic numbers and masses of the particles 
    are also given. The last argument is the gamma parameter.

    The function is not restricted to any fixed number of nodes in the integration.
    Any order of Gauss-Legendre quadrature can therefore be used.

    Units are assumed SI.

    """

    integral = 0
    theta = 0
    b_max = 10**(-9)
    E_com = conv.lab_to_com(E_lab, m1, m2)

    for i in range(len(weights)):

        # scale the b argument for the quadrature. 
        b = b_max/2*(abscissas[i] + 1)

        # find the rmin by solving g(r)^2 = 0
        r_min = fsolve(func=hf.g_squared, x0=10**(-10), args=(b, E_com, z1, z2)) 

        # calculate the inner integral, which is the scattering angle, using rmin as the lower bound
        for j in range(len(weights)):
            theta += weights[j] * hf.F((1/2*abscissas[j]+1/2), b, E_com, r_min, z1, z2)

        theta = 1/2*theta
        theta = np.pi - 4*b*theta

        # sum the stopping power terms
        integral += weights[i] * b*(np.sin(theta/2))**2

    integral = 2*np.pi*gamma*E_lab*b_max/2*integral

    return integral