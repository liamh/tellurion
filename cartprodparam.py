# Utilities to help with studies that repeat the same calculation with
# different combinations of input parameters, specifically, a
# Cartesian product of parameters. The function `callparamset()`
# generates the combinations, calls the function, collects the
# results in a data frame, and outputs to CSV.

import itertools
import io
import pandas as pd
from joblib import Parallel, delayed

####### Inclusive range

# In [29]: rangi(45)
# Out[29]: [45]

# In [30]: rangi(45,50)
# Out[30]: [45, 46, 47, 48, 49, 50]

# In [31]: rangi(45,50,5)
# Out[31]: [45, 50]

# In [32]: rangi(45,50,6)
# Out[32]: [45]
def rangi(start, stop=None, step=1):
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

####### Utilities for lists and dictionaries

# If argument is a list, return it; otherwise, make a singleton list
def ensurelist(arg):
    if type(arg) is list:
        return(arg)
    else:
        return([arg])

# Make a list of the dictionary values (no keys)
def dictvals(dict):
    return list(dict.values())

####### Outer product of lists with constraint

# Cartesian product of values from ranges, with optional constraint
# Ranges is a list of ranges, as arguments to rangi; singletons may be left as numbers
# cartprodrange([[250,325,25],[275],[45],[2021,2022],[1,12,3]])
# cartprodrange([[250,325,25],275,45,[2021,2022],[1,12,3]], lambda l: l[0] <= l[1])
def cartprodrange (ranges, constraint=None):
    return(cartprod(list(map(lambda range: rangi(*ensurelist(range)), ranges)), constraint))

# Cartesian products of sets
def cartprod (paramsets, constraint=None):
    ipr = itertools.product(*paramsets)
    if constraint==None:
        ret = [list(i) for i in ipr]
    else:
        ret = [list(i) for i in ipr if constraint(i)]
    return(ret)

####### Call function in parallel with combinations of arguments and collect results

# Call a function with integer arguments in every combination of
# specified ranges and create a CSV of the results

# Print to string https://stackoverflow.com/a/56103429/238405
def sprint(*args, end='', **kwargs):
    sio = io.StringIO()
    print(*args, **kwargs, end=end, file=sio)
    return sio.getvalue()

def nospaces(*objects):
    return sprint(*objects).replace(" ", "")

# Call the function with all argument combinations and collect the
# result in a CSV which is named in the return value. The function fn
# should return a dict with the same keys regardless of the arguments.
def callparamset(args, fn, constraint=None):
    arglists = cartprodrange(args, constraint)
    resultsll = Parallel(n_jobs=len(arglists))(delayed(fn)(*args) for args in arglists)
    df = pd.DataFrame(resultsll)
    filename = nospaces(fn.__name__,"-",args,".csv")
    df.to_csv(filename)
    return(filename)
