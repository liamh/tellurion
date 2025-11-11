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
from tellurion.core import nquant

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
    time: astropy.time.Time
    '''The date and time of the state'''
    aux: dict = dataclasses.field(default_factory=dict)
    '''Discrete attributes of the orbital state; these are attributes that have a finite set of discrete values'''
    _cartesian: u.Quantity = dataclasses.field(default=None, init=False)
    '''Internal storage for Cartesian coordinates'''
    _spherical: u.Quantity = dataclasses.field(default=None, init=False)
    '''Internal storage for spherical coordinates'''

    def __init__(self, time, aux=None, cartesian=None, spherical=None, pv=None):
        """Initialize PVT with either cartesian or spherical coordinates.

        Args:
            time: astropy.time.Time object
            aux: dictionary of auxiliary attributes
            cartesian: Cartesian state vector (position + velocity)
            spherical: Spherical state vector (r, theta, phi, vr, vtheta, vphi)
            pv: Legacy parameter name for cartesian (for backwards compatibility)
        """
        self.time = time
        self.aux = aux if aux is not None else {}

        # Handle legacy pv parameter
        if pv is not None and cartesian is None:
            cartesian = pv

        # Ensure only one coordinate system is provided
        if cartesian is not None and spherical is not None:
            raise ValueError("Cannot specify both cartesian and spherical coordinates")
        if cartesian is None and spherical is None:
            raise ValueError("Must specify either cartesian or spherical coordinates")

        self._cartesian = cartesian
        self._spherical = spherical

    @property
    def cartesian(self):
        """Lazy evaluation of Cartesian coordinates from spherical if needed."""
        if self._cartesian is None:
            self._cartesian = self._spherical_to_cartesian(self._spherical)
        return self._cartesian

    @cartesian.setter
    def cartesian(self, value):
        """Set Cartesian coordinates and clear cached spherical."""
        self._cartesian = value
        self._spherical = None

    @property
    def spherical(self):
        """Lazy evaluation of spherical coordinates from Cartesian if needed."""
        if self._spherical is None:
            self._spherical = self._cartesian_to_spherical(self._cartesian)
        return self._spherical

    @spherical.setter
    def spherical(self, value):
        """Set spherical coordinates and clear cached Cartesian."""
        self._spherical = value
        self._cartesian = None

    @property
    def pv(self):
        """Legacy property name for cartesian coordinates."""
        return self.cartesian

    @pv.setter
    def pv(self, value):
        """Legacy setter for cartesian coordinates."""
        self.cartesian = value

    def _cartesian_to_spherical(self, cartesian):
        """Convert Cartesian state vector to spherical coordinates.

        Args:
            cartesian: Cartesian state vector with position and velocity components

        Returns:
            spherical: Spherical state vector using sph() format
        """
        from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
        from astropy.coordinates import CartesianDifferential, SphericalDifferential

        # Extract position and velocity from the structured quantity
        pos = cartesian[_eph_pos]  # position 3-vector
        vel = cartesian[_eph_vel]  # velocity 3-vector

        # Create CartesianRepresentation with position
        cart_repr = CartesianRepresentation(x=pos[0], y=pos[1], z=pos[2])

        # Create CartesianDifferential with velocity
        cart_diff = CartesianDifferential(d_x=vel[0], d_y=vel[1], d_z=vel[2])

        # Add the differential to the representation
        cart_repr = cart_repr.with_differentials(cart_diff)

        # Convert to spherical representation (this handles both position and velocity)
        sph_repr = cart_repr.represent_as(SphericalRepresentation,
                                          differential_class=SphericalDifferential)

        # Extract the spherical differential
        sph_diff = sph_repr.differentials['s']  # 's' is the time unit key

        # Convert to the format expected by sph() function
        # SphericalRepresentation uses (lon, lat, distance) format
        # which corresponds to (right ascension, declination, distance)
        sphrepr = [sph_repr.lon, sph_repr.lat, sph_repr.distance]
        sphrate = [sph_diff.d_lon, sph_diff.d_lat, sph_diff.d_distance]

        # Use the existing sph() function to create the structured quantity
        return sph(sphrepr, sphrate, labels=['rtasc', 'decl', 'distance'],
                   unitlookup=astro.prefunits)

    def _spherical_to_cartesian(self, spherical):
        """Convert spherical coordinates to Cartesian state vector.

        Args:
            spherical: Spherical state vector in sph() format

        Returns:
            cartesian: Cartesian state vector with position and velocity components
        """
        from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
        from astropy.coordinates import CartesianDifferential, SphericalDifferential

        # Extract spherical position and velocity from the structured quantity
        # Assuming the sph() function creates fields: rtasc, decl, distance, rtasc_r, decl_r, distance_r
        rtasc = spherical['rtasc']          # right ascension (longitude)
        decl = spherical['decl']            # declination (latitude)
        distance = spherical['distance']    # radial distance
        rtasc_r = spherical['rtasc_r']      # d(right ascension)/dt
        decl_r = spherical['decl_r']        # d(declination)/dt
        distance_r = spherical['distance_r'] # d(distance)/dt

        # Create SphericalRepresentation with position
        # Note: SphericalRepresentation expects (lon, lat, distance)
        sph_repr = SphericalRepresentation(
            lon=rtasc,
            lat=decl,
            distance=distance
        )

        # Create SphericalDifferential with velocity
        sph_diff = SphericalDifferential(
            d_lon=rtasc_r,
            d_lat=decl_r,
            d_distance=distance_r
        )

        # Add the differential to the representation
        sph_repr = sph_repr.with_differentials(sph_diff)

        # Convert to Cartesian representation (this handles both position and velocity)
        cart_repr = sph_repr.represent_as(CartesianRepresentation,
                                          differential_class=CartesianDifferential)

        # Extract the Cartesian differential
        cart_diff = cart_repr.differentials['s']  # 's' is the time unit key

        # Extract position and velocity as 3-vectors
        position = [cart_repr.x, cart_repr.y, cart_repr.z]
        velocity = [cart_diff.d_x, cart_diff.d_y, cart_diff.d_z]

        # Use the existing pv() function to create the structured quantity
        return pv(position, velocity, unitlookup=astro.prefunits)

    def __getitem__(self, index):
        cart1d = util.ensure_1d(self.cartesian)
        t1d = util.ensure_1d(self.time)
        return PVT(time=t1d[index], cartesian=cart1d[index],
                   aux={k: v.split(' ')[index] for k, v in self.aux.items()})

    def __len__(self):
        if type(self.time) is np.ndarray:
            return len(self.time)
        else:
            return 1

    def copy(self):
        # Only copy the attribute that's currently defined to avoid unnecessary computation
        if self._cartesian is not None:
            return PVT(time=self.time.copy(), cartesian=self.cartesian.copy(), aux=self.aux.copy())
        else:
            return PVT(time=self.time.copy(), spherical=self.spherical.copy(), aux=self.aux.copy())

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
        self._cartesian = np.concatenate((util.ensure_1d(self.cartesian), util.ensure_1d(pvt.cartesian)))
        self._spherical = None  # Clear cached spherical coordinates
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
            ts = TimeSeries(time=tm, data=self.cartesian[_eph_pos], names=columnnames[1:])
        else:
            ts = TimeSeries(time=tm, data=self.cartesian, names=columnnames[1:])
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
            return np.hstack((self.cartesian.si[_eph_pos].value, self.cartesian.si[_eph_vel].value, \
                              self.time.to_array(time_format).reshape(1)))
        else:
            return np.hstack((self.cartesian.si[_eph_pos].value, self.cartesian.si[_eph_vel].value, \
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
    return type(obj) == PVT and ispv(obj.cartesian) and isdttm(obj.time)
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
        return nquant.structquant((position, velocity), (_eph_pos, _eph_vel), \
                             phystype=('length','speed'), unitlookup=unitlookup)

def sph(sphrepr, sphrate=None, labels=['rtasc','decl','distance'], unitlookup=astro.prefunits):
    '''Make a spherical coordinate set for position.
        Args:
         sphrepr = [ cylang: cylindrical angle (azimuth, longitude, right ascension),
                     sphang: spherical angle (elevation, latitude, declination),
                     distance: radial distance (range) ]
         sphrate = the time derivatives of sphrepr
        Returns:
         Structured quantity
    '''
    sphr = [u.Quantity(sphrepr[0]).to(unitlookup['angle']), \
            u.Quantity(sphrepr[1]).to(unitlookup['angle']), \
            u.Quantity(sphrepr[2]).to(unitlookup['length'])]
    if sphrate:
        labels_r = [sym+'_r' for sym in labels]
        return nquant.structquant(sphr+sphrate, labels + labels_r, \
                             phystype=('angle', 'angle', 'length', \
                                       'angular speed', 'angular speed', 'speed'), \
                            unitlookup=unitlookup)
    else:
        return nquant.structquant(sphr, labels, \
                             phystype=('angle', 'angle', 'length'), unitlookup=unitlookup)

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
        return PVT(time=obj[_eph_time], cartesian=pv(pos, vel))
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
            res = PVT(time=astro.abstime(0), cartesian=obj)
        else:
            # Add the specified time to the PV
            res = PVT(time=astro.abstime(item), cartesian=obj)
    elif ispvt(obj):
        if isdttm(item):
            # Replace the timestamp in the PVT
            res = PVT(time=item, cartesian=obj.cartesian)
        else:
            # Displace the timestamp in the PVT by the given relative time
            res = PVT(time=obj.time + item, cartesian=obj.cartesian)
    elif type(obj) is tuple:
        # Create a PVT from the three P, V, T
        res = PVT(time=obj[2], cartesian=pv(obj[0], obj[1]))
    elif type(obj) is list:
        if len(obj)>1:
            return obj[0].copy().concatenate(obj[1:])
        elif len(obj)==1:
            return obj[0]
    elif type(obj) is np.ndarray:
        return PVT(time=astro.abstime(astropy.time.Time(obj[6], format=astro.prefnumabstime).to_datetime()),
                   cartesian=pv(obj[0:3], obj[3:6], astro.siunits).to(astro.prefunits['posvel']))
    else:
        raise ValueError('Cannot make a PVT from this object')
    if aux: # only one attribute for now
        res[aux[0]] = aux[1]
    return res

def makepos(pos, unit=astro.prefunits['length']):
    '''Create a position vector or convert units'''
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

##################################################
#### Compare positions, velocities            ####
##################################################

def magdiff(a, b):
    '''Magnitude of the difference of two vectors'''
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])

def posdiff(a, b):
    return np.linalg.norm(a[_eph_pos]-b[_eph_pos])
