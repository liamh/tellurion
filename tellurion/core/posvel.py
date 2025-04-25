"""
Position, velocity and time sets in AstroPy
"""
import collections
import dataclasses
import datetime
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time
from astropy.timeseries import TimeSeries
import astropy.table.row
import astropy.coordinates as coord
from . import astro

##################################################
####   Constants used to define field names   ####
##################################################

# Provide attributes with default values https://stackoverflow.com/a/18348004/238405
# Maybe use dataclasses https://stackoverflow.com/q/47955263/238405
#PVT = collections.namedtuple('PVT', 'pv time aux')

@dataclasses.dataclass
class PVT:
    '''Orbital state vector as position (Cartesian 3-vector), velocity (Cartesian 3-vector), time, and a dictionary of discrete attributes; each field can have multiple rows, corresponding to an ephemeris'''
    pv: u.Quantity
    '''The orbital state vector as an astropy.units 6-vector with structured quantity of physical dimension length, speed'''
    time: astropy.time.Time
    '''The date and time of the state'''
    aux: dict = dataclasses.field(default_factory=dict)
    '''Discrete attributes of the orbital state; these are attributes that have a finite set of discrete values'''

    def copy(self):
        return PVT(self.pv.copy(), self.time.copy(), self.aux.copy())

    def vcat(self, pvt):
        if self.pv.shape==():
            reshapepv1 = np.reshape(self.pv, (1,))
        else:
            reshapepv1 = self.pv
        if pvt.pv.shape==():
            reshapepv2 = np.reshape(pvt.pv, (1,))
        else:
            reshapepv2 = pvt.pv
        self.pv = np.concatenate((reshapepv1, reshapepv2))

        if type(self.time.value) is list:
            t1 = self.time.value
        else:
            t1 = [self.time.value]
        if type(pvt.time.value) is list:
            t2 = pvt.time.value
        else:
            t2 = [pvt.time.value] # need here pvt.time.to_string()
        breakpoint() # times are as datetimes
        self.time = astropy.time.Time(t1 + t2)
        dl = [self.aux, pvt.aux]
        self.aux ={k: [d[k] for d in dl] for k in dl[0]}
        return self

_eph_time = 'time'
_eph_pos = 'position'
_eph_vel = 'velocity'
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]
_eph_pos_xyz = ('px','py','pz')

# These should be conditional on the units used
_pos_format = '10.3f'
_vel_format = '10.6f'

def tsephem(data, times=None, columnnames=_ephemeris_columns, \
            pvformats=(_pos_format, _vel_format)):
    '''Create an ephemeris table (timeseries) from the posvel and times or PVT data'''
    if times==None:
        (data, times, auxes) = zip(*data)
    datdict = {columnnames[1]: [d[_eph_pos] for d in data], \
               columnnames[2]: [d[_eph_vel] for d in data]}
    ts = TimeSeries(time=times, data=datdict)
    ts[columnnames[1]].info.format = pvformats[0]
    ts[columnnames[2]].info.format = pvformats[1]
    return ts

def posxyz(ephem):
    '''Position only, separate x, y, z, components'''
    txyz = TimeSeries(time=ephem['time'], data=ephem[_eph_pos])
    txyz.rename_columns(('col0','col1','col2'), _eph_pos_xyz)
    return txyz

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

# A PVT consists of a tuple a posvel (as defined by ispv()) and an astropy.time.Time
def pvt(obj, item=None, aux=None):
    '''Return a tuple of posvel and time from a variety of sources.
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
            res = PVT(obj[0], item)
        else:
            # Displace the timestamp in the PVT by the given relative time
            res = PVT(obj[0], obj[1] + item)
    elif type(obj) is tuple:
        # Create a PVT from the three P, V, T
        res = PVT(pv(obj[0], obj[1]), obj[2])
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

def pvtsijd(pvt, units=astro.prefunits):
    '''Create an array of length 7 with position and velocity in SI units and Julian date, or create a PVT from an array of length 7'''
    if type(pvt) is PVT:
        return np.concatenate((pvt.pv.si[_eph_pos].value, pvt.pv.si[_eph_vel].value, \
                               np.array([pvt.time.to_value('jd')])))
    elif type(pvt) is np.ndarray:
        return PVT(pv(pvt[0:3], pvt[3:6], astro.siunits).to(units['posvel']), \
                   astro.abstime(astropy.time.Time(pvt[6], format='jd').to_datetime()))

##################################################
#### Compare positions, velocities            ####
##################################################

def magdiff(a, b):
    '''Magnitude of the difference of two vectors'''
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])

def posdiff(a, b):
    return np.linalg.norm(a[_eph_pos]-b[_eph_pos])
