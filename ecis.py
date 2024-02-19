# ECIS exploratory computation in stages
from munch import *
class Ecis(Munch):
    pass

def newtree(name=None, value=None):
    ret = Ecis()
    if name is not None:
        ret.update({name: value})
    return(ret)
