"""
Position, velocity and time sets in AstroPy
"""
import collections
import dataclasses
import datetime
import funcy
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time
from astropy.timeseries import TimeSeries
import astropy.table.row
import astropy.coordinates as coord
from tellurion.core import util
from tellurion.core import astro

##################################################
####   Constants used to define field names   ####
##################################################

_eph_time = 'time'
_eph_pos = 'position'
_eph_vel = 'velocity'
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]
_ephemeris_columns_pos_xyz = [_eph_time, 'px','py','pz']

# These should be conditional on the units used
_pos_format = '10.3f'
_vel_format = '10.6f'

# Provide attributes with default values https://stackoverflow.com/a/18348004/238405
# Maybe use dataclasses https://stackoverflow.com/q/47955263/238405
#PVT = collections.namedtuple('PVT', 'pv time aux')

@dataclasses.dataclass
class PVT(collections.abc.Sequence):
    '''Orbital state vector as position (Cartesian 3-vector), velocity (Cartesian 3-vector), time, and a dictionary of discrete attributes; each field can have multiple rows, corresponding to an ephemeris'''
    pv: u.Quantity
    '''The orbital state vector as an astropy.units 6-vector with structured quantity of physical dimension length, speed'''
    time: astropy.time.Time
    '''The date and time of the state'''
    aux: dict = dataclasses.field(default_factory=dict)
    '''Discrete attributes of the orbital state; these are attributes that have a finite set of discrete values'''

    def __getitem__(self, index):
        pv1d = util.ensure_1d(self.pv)
        t1d = util.ensure_1d(self.time)
        return PVT(pv1d[index], t1d[index],
                   {k: v.split(' ')[index] for k, v in self.aux.items()})

    def __len__(self):
        if type(self.time) is np.ndarray:
            return len(self.time)
        else:
            return 1

    def copy(self):
        return PVT(self.pv.copy(), self.time.copy(), self.aux.copy())

    def timeorder(self):
        '''Sort the PVT in increasing time order'''
        return pvt(sorted(self, key=lambda x: x.time))

    def concatenate(self, pvt):
        '''Concatenate rows from another PVT onto the end of this PVT'''
        if type(pvt) is list:
            if pvt:
                return self.concatenate(pvt[0]).concatenate(pvt[1:])
            else:
                return self
        self.pv = np.concatenate((util.ensure_1d(self.pv), util.ensure_1d(pvt.pv)))
        self.time = astro.abstime([self.time, pvt.time])
        self.aux = funcy.merge_with(' '.join, self.aux, pvt.aux)
        return self

    def merge(self, pvt):
        '''Merge the two PVTs and put in time order'''
        return self.copy().concatenate(pvt).timeorder()

    def ephemeris(self, elapsed=True, reftime='epoch', \
                  columnnames=_ephemeris_columns, pvformats=(_pos_format, _vel_format)):
        '''An AstroPy time series of PVT with `aux` variables (if
        any). If `elapsed=True` (default), then add a column that
        gives the elapsed time from the previous time step.'''
        if self.aux:
            columnnames = _ephemeris_columns_pos_xyz
        if self.time.isscalar:
            tm = astropy.time.Time([self.time.to_string()])
        else:
            tm = self.time
        if len(columnnames)==4:
            ts = TimeSeries(time=tm, data=self.pv[_eph_pos], names=columnnames[1:])
        else:
            ts = TimeSeries(time=tm, data=self.pv, names=columnnames[1:])
        for key in self.aux:
            ts[key] = self.aux[key].split(' ')
        if len(columnnames)==4:
            ts[columnnames[1]].info.format = pvformats[0]
            ts[columnnames[2]].info.format = pvformats[0]
            ts[columnnames[3]].info.format = pvformats[0]
        else:
            ts[columnnames[1]].info.format = pvformats[0]
            ts[columnnames[2]].info.format = pvformats[1]
        if elapsed:
            elapsed = [dt.quantity_str for dt in np.diff(tm)]
            elapsed.insert(0,'')
            ts.add_column(elapsed, index=1, name='elapsed')
        if reftime is not None:
            astro.fromtime(ts, reftime=reftime, copy=False)
        return ts

    def to_array(self, time_format=astro.prefnumabstime):
        '''Convert the PVT into a 7-column np.ndarray using SI units and the preferred time format ('mjd' default); `aux` values are not included'''
        if self.time.shape==():
            return np.hstack((self.pv.si[_eph_pos].value, self.pv.si[_eph_vel].value, \
                              self.time.to_array(time_format).reshape(1)))
        else:
            return np.hstack((self.pv.si[_eph_pos].value, self.pv.si[_eph_vel].value, \
                              self.time.to_array(time_format).reshape(-1,1)))

##################################################
####   Tests for posvel and related types     ####
##################################################

def isq3vec(obj, physdim):
    '''Is a 3-vector u.Quantity with the specified physical dimension'''
    return type(obj) is u.Quantity \
        and u.get_physical_type(obj) == u.get_physical_type(physdim) \
        and np.size(obj)==3

def ispv(obj):
    return type(obj) == u.Quantity and obj.dtype.names is not None \
        and _eph_pos in obj.dtype.names and _eph_vel in obj.dtype.names \
        and isq3vec(obj[_eph_pos],'length') and isq3vec(obj[_eph_vel],'speed')

def isdttm(obj):
    return type(obj) is astropy.time.Time

def isreltime(obj):
    """Object is a relative time: is a u.Quantity with physical type 'time'"""
    return type(obj) is u.Quantity and u.get_physical_type(obj) == u.get_physical_type('time')

def ispvt(obj):
    return type(obj) == PVT and ispv(obj.pv) and isdttm(obj.time)
#or (type(obj) == tuple and len(obj) == 2 \
#        and ispv(obj[0]) and isdttm(obj[1]))

def isephem(ts):
    '''Argument is an ephemeris table'''
    return type(ts) is TimeSeries \
        and all([k in ts.keys() for k in _ephemeris_columns])

def isephrow(row):
    '''Argument is a row of an ephemeris table'''
    return type(row) is astropy.table.row.Row \
        and all([row.keys().__contains__(k) for k in _ephemeris_columns])

def ispvter(obj):
    return ispvt(obj) or isephrow(obj)

##################################################
####   Make posvel and related types          ####
##################################################

def pv(position, velocity, unitlookup=astro.prefunits):
    '''Make a posvel from separate position and velocity; if argument `position` is a posvel, then convert units'''
    if ispv(position):
        return(position.to(unitlookup['length']))
    else:
        return astro.makesq((position, velocity), (_eph_pos, _eph_vel), \
                            phystype=('length','speed'), unitlookup=unitlookup)

##################################################
#### Dates, times, and PVT                    ####
##################################################

# To add timezone to datetime
#import pytz
#def utcdt(datetime):
#    return pytz.utc.localize(datetime)

def pvt(obj, item=None, aux=None):
    '''Return an instance of PVT (position, velocity and time), from a
    variety of sources. Sequence in time can be length 1 or more.

    '''
    if isephrow(obj):
        # PVT from an ephemeris row
        pos = obj[_eph_pos]
        vel = obj[_eph_vel]
        res = PVT(pv(pos, vel), obj[_eph_time])
    elif type(obj) is TimeSeries:
        # Select a row from an ephemeris by index, absolute time, or relative time
        if item == None:   # Return the last row
            item = -1
        try:
            row = obj[item]
        except:
            row = obj.loc[astro.abstime(obj[0]['time'])]
        return pvt(row, aux = aux)
    elif ispv(obj):
        if item==None:
            # Add the current time to the PV
            res = PVT(obj, astro.abstime(0))
        else:
            # Add the specified time to the PV
            res = PVT(obj, astro.abstime(item))
    elif ispvt(obj):
        if isdttm(item):
            # Replace the timestamp in the PVT
            res = PVT(obj.pv, item)
        else:
            # Displace the timestamp in the PVT by the given relative time
            res = PVT(obj.pv, obj.time + item)
    elif type(obj) is tuple:
        # Create a PVT from the three P, V, T
        res = PVT(pv(obj[0], obj[1]), obj[2])
    elif type(obj) is list:
        if len(obj)>1:
            return obj[0].copy().concatenate(obj[1:])
        elif len(obj)==1:
            return obj[0]
    elif type(obj) is np.ndarray:
        return PVT(pv(obj[0:3], obj[3:6], astro.siunits).to(astro.prefunits['posvel']), \
                   astro.abstime(astropy.time.Time(obj[6], format=astro.prefnumabstime).to_datetime()))
    else:
        raise ValueError('Cannot make a PVT from this object')
    if aux: # only one attribute for now
        res[aux[0]] = aux[1]
    return res

def makepos(pos, unit=astro.prefunits['length']):
    '''Create a position vector or convert units'''
    if type(pos) is coord.SkyCoord:
        return makepos(pos.cartesian, unit)
    if type(pos) is coord.representation.cartesian.CartesianRepresentation:
        return makepos(pos.xyz, unit)
    if type(pos) is tuple:
        return u.Quantity(pos, unit)
    if u.get_physical_type(unit) == 'length':
        if isq3vec(pos, 'length'):
            return pos.to(unit)
        elif type(pos) is u.Quantity:
            raise ValueError('Argument does not represent a position 3-vector')
        else:
            return pos*unit
    else:
        raise ValueError('Unit does not represent a position')

def makept(pos, dttm, frame='gcrs'):
    '''Make a SkyCoord GCRS location'''
    return coord.SkyCoord(coord.CartesianRepresentation(pos), obstime=dttm, frame=frame)

##################################################
#### Compare positions, velocities            ####
##################################################

def magdiff(a, b):
    '''Magnitude of the difference of two vectors'''
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])

def posdiff(a, b):
    return np.linalg.norm(a[_eph_pos]-b[_eph_pos])
