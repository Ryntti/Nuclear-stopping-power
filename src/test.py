import integral
import numpy as np
import constants as c
import conversions as conv
import semiempirical as sem
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss 

# order of the Gauss-Legendre approximationn. number of abscissas
n = 8
# get the nodes and weights
abscissas, weights = leggauss(n)

# get the initial projectile energies in center of mass coordinates
E_com = [conv.lab_to_com(x, c.m1[0], c.m2[0]) for x in c.E_lab]

inter = integral.sp(abscissas, weights, c.E_lab[0], c.z1[0], c.z2[0], c.m1[0], c.m2[0], c.gamma[0])
print(conv.J_to_eV([inter]))










