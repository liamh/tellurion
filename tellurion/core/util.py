import numpy as np
import astropy.units as u

##################################################
#### Compare vectors by their magnitudes      ####
##################################################

def magdiff(a, b):
    """Magnitude of the difference of two vectors"""
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])
