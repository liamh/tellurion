import astropy.units as u
import numpy as np

##################################################
#### Compare vectors by their magnitudes      ####
##################################################

def magdiff(a, b):
    """Magnitude of the difference of two vectors"""
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])
