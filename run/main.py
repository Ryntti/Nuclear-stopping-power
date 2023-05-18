import sys
sys.path.insert(0, "C:\\Studies\\NMSC\\piirila_reino_nmsc_project1\\src")

import integral
import constants as c
import conversions as conv
import semiempirical as sem
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss 

if __name__ == "__main__":

    # order of the Gauss-Legendre approximationn. number of abscissas
    n = 8
    # get the nodes and weights
    abscissas, weights = leggauss(n)
    
    # get the initial kinetic energies of the system in center of mass coordinates
    E_com = [conv.lab_to_com(x, c.m1[0], c.m2[0]) for x in c.E_lab]

    stopping_power1 = []
    semesp1 = []
    for i in range(len(E_com)):
        stopping_power1.append(integral.sp(abscissas, weights, c.E_lab[i], c.z1[0], c.z2[0], c.m1[0], c.m2[0], c.gamma[0]))
        semesp1.append(sem.big_sn(c.E_lab_eV[i]*10**(-3),c.m1[0], c.m2[0], c.z1[0], c.z2[0]))

    stopping_power2 = []
    semesp2 = []
    for i in range(len(E_com)):
        stopping_power2.append(integral.sp(abscissas, weights, c.E_lab[i], c.z1[1], c.z2[1], c.m1[1], c.m2[1], c.gamma[1]))
        semesp2.append(sem.big_sn(c.E_lab_eV[i]*10**(-3), c.m1[1], c.m2[1], c.z1[1], c.z2[1]))

    plt.figure()
    stopping_power1 = conv.J_to_eV(stopping_power1)
    plt.title('Stopping power as a function of initial projectile energy in lab coordinates')
    plt.loglog(c.E_lab_eV, stopping_power1, 'r-', label='stopping power')
    plt.loglog(c.E_lab_eV, semesp1, 'b-', label='semiempirical stopping power')
    plt.xlabel('Initial projectile kinetic energy (eV)')
    plt.ylabel('Nuclear stopping power (eV*m^2)')
    plt.legend()
    plt.savefig('Stopping_power1.png')
    plt.show()

