import math
import numpy as np
import conversions as conv
# A module for constant and other parameter declarations 

# constants:
e0 = 8.854187813 * 10**(-12)                      # vacuum permittivity (SI)
q =  1.602176634 * 10**(-19)                      # elementary charge (SI)
proton_mass = 1.672621924 * 10**(-27)             # kg (SI)
Si28_mass = 4.645681117 * 10**(-26)               # kg (SI)
Au197_mass = 3.270709053 * 10**(-25)              # kg (SI)

# parameters:
alpha = [0.1818, 0.5099, 0.2802, 0.02817]         # alpha parameters (dimensionless)
beta = [3.2, 0.9423, 0.4028, 0.2016]              # beta parameters (dimensionless)
z1 = [1, 28]                                      # Atomic numbers of the projectiles (1st hydrogen, 2nd Si)
z2 = [28, 197]                                    # Atomic numbers of the targets (1st Si, 2nd Au)
m1 = [proton_mass, Si28_mass]                     # Masses of the projectiles (1st proton, 2nd Si)
m2 = [Si28_mass, Au197_mass]                      # Masses of the targets (1st Si, 2nd Au)

# gamma parameters:
gamma = [4*m1[0]*m2[0]/(m1[0]+m2[0])**2, 4*m1[1]*m2[1]/(m1[1]+m2[1])**2]

# initial projectile kinetic energy values:
power = math.log10(5*10**6)
E_lab_eV = np.logspace(1, power, 20)              # Initial lab coordinate energies of the projectiles (eV)
E_lab = conv.eV_to_J(E_lab_eV)                    # Same energies in Joules
