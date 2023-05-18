import math 
import conversions as conv

def eps(E_lab, m1, m2, z1, z2):
    """
    A function for the epsilon parameter featured in the semiempirical 
    S_n formula. 
    Units of E_lab are assumed keV.
    """
    return 32.53*m2*E_lab/(z1*z2*(m1+m2)*(z1**(0.23) + z2**(0.23)))

def sn(E_lab, m1, m2, z1, z2):
    """
    A helper function for the big_sn function. 
    Units of E_lab assumed keV.

    """
    ep = eps(E_lab, m1, m2, z1, z2)
    if ep <= 30:
        return math.log(1 + 1.138*ep) / (2*(ep +  0.01321*ep**0.21226 + 0.19593*ep**0.5))
    elif ep > 30:
        return math.log(ep)/(2*ep)
    
def big_sn(E_lab: float, m1: float, m2: float, z1: int, z2: int):
    """
    The semiempirical formula of stopping power as a function. 
    Units of E_lab assumed keV.

    """
    # E_lab = conv.J_to_keV(E_lab)
    return 8.462*10**(-15)*z1*z2*m1 / ((m1+m2)*(z1**0.23 + z2**0.23))*sn(E_lab, m1, m2, z1, z2)







