# Python utilities

import itertools

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
# Ranges is a list of ranges, as arguments to rangi
# coutprod([[250,325,25],[275],[45],[2021,2022],[1,12,3]])
# coutprod([[250,325,25],[275],[45],[2021,2022],[1,12,3]], lambda l: l[0] <= l[1])
def coutprod (ranges, constraint=None):
    # iparg = [rangi(*ranges[0]), rangi(*ranges[1]), rangi(*ranges[2]), rangi(*ranges[3]), rangi(*ranges[4])]
    iparg = list(map(lambda range: rangi(*range), ranges))
    ipr = itertools.product(*iparg)
    if constraint==None:
        ret = [list(i) for i in ipr]
    else:
        ret = [list(i) for i in ipr if constraint(i)]
    return(ret)
