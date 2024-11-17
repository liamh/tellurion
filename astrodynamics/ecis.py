""" ECIS exploratory computation in stages

Intermediate values of a computation are saved to a
structure/dict. New computations add to previously generated results,
possibly with a new "branch" from the tree.
"""

import munch
import warnings
import collections.abc

ecisdefault = 'ecdefault'

class Ecis(munch.Munch):
    def __init__(self, defaultdict=None):
        if defaultdict is not None:
            self.update({ecisdefault: defaultdict})

def newtree(name=None, value=None, defaultdict=None):
    """
    Create a new ECIS tree for computations.
    """
    ret = Ecis(defaultdict)
    if name is not None:
        ret.update({name: value})
    return(ret)

def _getobjs(tree, classtype):
    return([key for key in tree if isinstance(tree[key], classtype)])

def thingofclass(object, classtype):
    """
    Return the item in the tree if it is of the requested
    classtype.  If there is more than one, the last one in the tree is
    returned.  Performs a recursive search in branches (ecis
    instances).  Returns a list of [item, dict].
    """
    if issubclass(type(object), classtype):
        thing = [object, None]
    elif isinstance(object, collections.abc.Iterable) and all(issubclass(type(ele), classtype) for ele in object):
        thing = [[thingofclass(ele,classtype)[0] for ele in object], None]
    elif type(object) is Ecis:
        objs = _getobjs(object, classtype)
        match len(objs):
           case 0:
              subtrees = _getobjs(object, Ecis)
              if len(subtrees)==0:
                  thing = [None, object]
                  #raise ValueError(f"No `{classtype.__name__}` objects found in tree")
              else:
                  for st in subtrees:
                      branch = object[st]
                      thing = [thingofclass(branch, classtype)[0], branch]
           case 1:
              thing = [object[objs[0]], object]
           case _:
              thing = [object[objs[-1]], object]
              warnings.warn(f"Multiple `{classtype.__name__}` objects defined in tree, using `{objs[-1]}`")
    else:
        thing=None
        raise ValueError(f"Not a `{classtype.__name__}` or tree")
    return(thing)


## Would be nice to have a "collapse" function to eliminate
## intermediate values that are not needed again and have no branches.
