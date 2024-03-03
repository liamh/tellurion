# ECIS exploratory computation in stages
from munch import *
import warnings

class Ecis(Munch):
    pass

def newtree(name=None, value=None):
    ret = Ecis()
    if name is not None:
        ret.update({name: value})
    return(ret)

def getobjs(tree, classtype):
    return([key for key in tree if isinstance(tree[key], classtype)])

def thingofclass(tree, classtype):
    if type(tree) is classtype:
        ret = tree
    elif type(tree) is Ecis:
        objs = getobjs(tree, classtype)
        match len(objs):
           case 0:
              subtrees = getobjs(tree, Ecis)
              if len(subtrees)==0:
                  thing = [None, tree]
                  #raise ValueError(f"No `{classtype.__name__}` objects found in tree")
              else:
                  for st in subtrees:
                      branch = tree[st]
                      thing = [thingofclass(branch, classtype)[0], branch]
           case 1:
              thing = [tree[objs[0]], tree]
           case _:
              thing = [tree[objs[-1]], tree]
              warnings.warn(f"Multiple `{classtype.__name__}` objects defined in tree, using `{objs[-1]}`")
    else:
        thing=Nothing
        raise ValueError(f"Not a `{classtype.__name__}` or tree")
    return(thing)


## Would be nice to have a "collapse" function to eliminate
## intermediate values that are not needed again and have no branches.
