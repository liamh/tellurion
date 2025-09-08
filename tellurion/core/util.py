import collections
import dataclasses
import datetime
import funcy
import numpy as np
import astropy.units as u
import astropy.time
from tellurion.core import astro

def ensure_1d(obj):
    '''Make object with shape `(n,)` where n=1 if obj is a scalar'''
    if obj.shape==():
        return obj.reshape(1,)
    else:
        return obj

class QuantT(collections.abc.Sequence):
    '''Base class for quantities and time, like position-velocity-time or observations and time.'''
    q: u.Quantity
    time: astropy.time.Time

    def __getitem__(self, index):
        q1d = ensure_1d(self.q)
        t1d = ensure_1d(self.time)
        return QuanT(q1d[index], t1d[index])

    def __len__(self):
        if type(self.time) is np.ndarray:
            return len(self.time)
        else:
            return 1

    def copy(self):
        return QuantT(self.q.copy(), self.time.copy())

    def timeorder(self):
        '''Sort in increasing time order'''
        return QuantT(sorted(self, key=lambda x: x.time))

    def concatenate(self, qt):
        '''Concatenate rows from another QT onto the end of this QT'''
        if type(qt) is list:
            if qt:
                return self.concatenate(qt[0]).concatenate(qt[1:])
            else:
                return self
        self.q = np.concatenate((ensure_1d(self.q), ensure_1d(qt.q)))
        self.time = astro.abstime([self.time, qt.time])
        return self

    def merge(self, qt):
        '''Merge the two QTs and put in time order'''
        return self.copy().concatenate(qt).timeorder()
