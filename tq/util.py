import numpy as np
import itertools

####### List and dict utilities

def ensurelist(arg):
    '''Ensure that the argument is a list; if it's not, wrap it in a list, otherise, return it.'''
    if arg is None:
        return([])
    if type(arg) is list:
        return(arg)
    return([arg])

def dictvals(dict):
    '''Dictionary values as a list.'''
    return list(dict.values())

def listnpa(thing):
    '''The argument is a list or numpy ndarray.'''
    return type(thing) is list or type(thing) is np.ndarray

####### Inclusive range

# Inclusive range
# In [29]: rangi(45)
# Out[29]: [45]

# In [30]: rangi(45,50)
# Out[30]: [45, 46, 47, 48, 49, 50]

# In [31]: rangi(45,50,5)
# Out[31]: [45, 50]

# In [32]: rangi(45,50,6)
# Out[32]: [45]
def rangi(start, stop=None, step=1):
    '''Inclusive range for numeric values.'''
    if (stop==None): # Return a singleton
        return ([start])
    # Find semi-exclusive range
    rg = range(start, stop, step)
    # Add the endpoint if it is not after the range
    if rg[-1]+step == stop:
        rgs = itertools.chain(rg, (stop,))
    else:
        rgs = rg
    return ([i for i in rgs])
