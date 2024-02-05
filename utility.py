# Utilities for calling Python functions with combinations of arguments

import io
import itertools
from joblib import Parallel, delayed
import pandas as pd


####### List and dict utilities

def ensurelist(arg):
    if type(arg) is list:
        return(arg)
    else:
        return([arg])

def dictvals(dict):
    return list(dict.values())

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

####### Outer product of lists with constraint

# Outer (Cartesian) product of values from ranges, with optional constraint
#  ranges:     a list of ranges, as arguments to rangi; singletons may be left as numbers
#  constraint: `None` or a boolean function of one argument, a list of values;
#              True indicates it is included, False that it is not
# Examples
#  coutprod([[250,325,25],[275],[45],[2021,2022],[1,12,3]])
#  coutprod([[250,325,25],275,45,[2021,2022],[1,12,3]], lambda l: l[0] <= l[1])
def coutprod (ranges, constraint=None):
    iparg = list(map(lambda range: rangi(*ensurelist(range)), ranges))
    ipr = itertools.product(*iparg)
    if constraint==None:
        ret = [list(i) for i in ipr]
    else:
        ret = [list(i) for i in ipr if constraint(i)]
    return(ret)

####### Call a function in parallel with integer arguments in every combination of specified ranges

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
def allargcomb(args,constraint,fn):
    arglists = coutprod(args, constraint)
    resultsll = Parallel(n_jobs=len(arglists))(delayed(fn)(*args) for args in arglists)
    df = pd.DataFrame(resultsll)
    filename = nospaces(fn.__name__,"-",args,".csv")
    df.to_csv(filename)
    return(filename)
