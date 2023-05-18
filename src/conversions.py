# A function which transforms the units of a list of energy values from eV to Joules
def eV_to_J(x: list):
    conv_factor = 1.602176634*10**(-19)
    return [t*conv_factor for t in x]

def keV_to_J(x: list):
    conv_factor = 1.602176634*10**(-16)
    return [t*conv_factor for t in x]


def J_to_eV(x: list):
    conv_factor = 6.241509074*10**18
    return [t*conv_factor for t in x]

def J_to_keV(x: float):
    conv_factor = 6.241509074*10**15
    return x*conv_factor

# A function which transforms the energy of the projectile in a binary collision from 
# laboratory coordinates to center of mass coordinates 
def lab_to_com(x, m1, m2):
    return x*m2/(m1+m2)







