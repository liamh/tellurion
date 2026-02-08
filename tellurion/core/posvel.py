"""
Position, velocity and time sets in AstroPy
"""
import abc
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
from tellurion.astro import units
from tellurion.astro import quantity_utils as quant
import tellurion.astro.time as atime
from tellurion.core import pvhelper

# Make the scalar structured quantity a singleton vector
u.Quantity.tovector = lambda self: u.Quantity([self]) if self.isscalar else self

##################################################
####  Base class for position-bearing objects ####
##################################################

@dataclasses.dataclass
class PositionBase(abc.ABC):
    """Base class for objects with cartesian and/or spherical position information.

    Provides lazy conversion between Cartesian and spherical coordinates,
    with the actual conversion logic implemented by subclasses depending
    on whether velocity information is present.
    """
    time: astropy.time.Time = None
    """The date and time of the position (optional)"""
    aux: dict = dataclasses.field(default_factory=dict)
    """Discrete attributes (labels, flags, metadata, etc.)"""
    _cartesian: u.Quantity = dataclasses.field(default=None, init=False, repr=False)
    """Internal storage for Cartesian coordinates"""
    _spherical: u.Quantity = dataclasses.field(default=None, init=False, repr=False)
    """Internal storage for spherical coordinates"""

    def __init__(self, time, aux=None, cartesian=None, spherical=None, pv=None):
        """Initialize PositionVelocityT with either cartesian or spherical coordinates.

        Args:
            time: astropy.time.Time object
            aux: dictionary of auxiliary attributes
            cartesian: Cartesian state vector (position + velocity) or array of state vectors
            spherical: Spherical state vector (r, theta, phi, vr, vtheta, vphi)
            pv: Legacy parameter name for cartesian (for backwards compatibility)
        """
        self.time = time
        self.isscalar = self.time.isscalar
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

    def __len__(self):
        if self.isscalar:
            raise TypeError(f"object of type '{type(self).__name__}' has no len()")
        else:
            return len(self.time)

    def __getitem__(self, index):
        """Get a single time step or slice."""
        if self.isscalar:
            raise TypeError(f"'{type(self).__name__}' scalar object is not subscriptable")
        else:
            return type(self)( # or equivalently self.__class__
                time=self.time[index],
                cartesian=self.cartesian[index],
                aux={k: v.split(' ')[index] for k, v in self.aux.items()} if self.aux else {}
            )

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
    def position_vector(self):
        """Get just the position 3-vector, regardless of whether velocity is present.

        Returns:
            u.Quantity: 3-vector with position components
        """
        cart = self.cartesian
        # Check if it's a structured quantity with position/velocity fields
        if hasattr(cart, 'dtype') and cart.dtype.names and pvhelper._eph_pos in cart.dtype.names:
            return cart[pvhelper._eph_pos]
        else:
            return cart

    @property
    def has_velocity(self):
        """Check if velocity information is present.

        Returns:
            bool: True if velocity data is available
        """
        cart = self.cartesian
        return (hasattr(cart, 'dtype') and cart.dtype.names and
                pvhelper._eph_vel in cart.dtype.names)

    @abc.abstractmethod
    def _cartesian_to_spherical(self, cartesian):
        """Convert Cartesian to spherical coordinates.

        Implemented by subclasses to handle position-only vs position+velocity.

        Args:
            cartesian: Cartesian coordinates (format depends on subclass)

        Returns:
            Spherical coordinates in sph() format
        """
        pass

    @abc.abstractmethod
    def _spherical_to_cartesian(self, spherical):
        """Convert spherical to Cartesian coordinates.

        Implemented by subclasses to handle position-only vs position+velocity.

        Args:
            spherical: Spherical coordinates in sph() format

        Returns:
            Cartesian coordinates (format depends on subclass)
        """
        pass

    @abc.abstractmethod
    def copy(self):
        """Create a copy of this object.

        Implemented by subclasses to return the correct type.
        """
        pass

    @abc.abstractmethod
    def _make_instance(self, time, cartesian, aux):
        """Factory method to create a new instance of the correct type.

        Args:
            time: astropy.time.Time object
            cartesian: Cartesian coordinates
            aux: Auxiliary attributes dictionary

        Returns:
            New instance of the appropriate subclass
        """
        pass

    def timeorder(self):
        """Sort in increasing time order.

        Returns:
            New instance of same type, sorted by time
        """
        if self.time is None or self.time.isscalar:
            return self.copy()
        sorted_items = sorted(self, key=lambda x: x.time)
        return self._from_sorted_list(sorted_items)

    def _from_sorted_list(self, sorted_list):
        """Create instance from a sorted list of single-time instances."""
        if not sorted_list:
            raise ValueError("Cannot create from empty list")
        if len(sorted_list) == 1:
            return sorted_list[0]

        # Concatenate all items
        result = sorted_list[0].copy()
        for item in sorted_list[1:]:
            result = result.concatenate(item)
        return result

    def concatenate(self, other):
        """Concatenate rows from another object onto the end of this object.

        Args:
            other: Another instance of the same type, or a list of instances

        Returns:
            self (modified in place)
        """
        if type(other) is list:
            if other:
                return self.concatenate(other[0]).concatenate(other[1:])
            else:
                return self

        # Concatenate cartesian coordinates
        self._cartesian = np.concatenate((self.cartesian.tovector(),
                                          other.cartesian.tovector()))
        self._spherical = None  # Clear cached spherical coordinates

        # Concatenate time
        if self.time is None and other.time is None:
            self.time = None
        elif self.time is None:
            self.time = other.time
        elif other.time is None:
            pass  # Keep self.time
        else:
            self.time = atime.abstime([self.time, other.time])

        self.aux = funcy.merge_with(' '.join, self.aux, other.aux) # Merge aux attributes
        self.isscalar = False
        return self

    def merge(self, other):
        """Merge with another object and put in time order.

        Args:
            other: Another instance of the same type

        Returns:
            New instance of same type, merged and sorted
        """
        return self.copy().concatenate(other).timeorder()

    def to_array(self, time_format=atime.prefnumabstime):
        """Convert to a numpy array using SI units.

        Args:
            time_format: Format for time column (default: 'mjd')

        Returns:
            np.ndarray: For PositionT: [px, py, pz, time] (4 columns)
                       For PositionVelocityT: [px, py, pz, vx, vy, vz, time] (7 columns)

        Note: aux attributes are not included in the array output
        """
        pos = self.position_vector

        # Check if we have velocity
        has_velocity = self.has_velocity

        if self.time is None:
            # No time information
            if has_velocity:
                vel = self.cartesian[pvhelper._eph_vel]
                return np.hstack((pos.si.value, vel.si.value))
            else:
                return pos.si.value
        else:
            # Has time information
            if self.time.shape == ():
                time_array = self.time.to_array(time_format).reshape(1)
            else:
                time_array = self.time.to_array(time_format).reshape(-1, 1)

            if has_velocity:
                vel = self.cartesian[pvhelper._eph_vel]
                return np.hstack((pos.si.value, vel.si.value, time_array))
            else:
                return np.hstack((pos.si.value, time_array))

    def ephemeris(self, elapsed=True, reftime='epoch', columnnames=None,
                  pvformats=(pvhelper._pos_format, pvhelper._vel_format)):
        """Create an AstroPy time series.

        For PositionT: Creates a time series with position only.
        For PositionVelocityT: Creates a time series with position and velocity.
        Both include any `aux` attributes.

        Args:
            elapsed: If True (default), add a column with elapsed time from previous step
            reftime: Reference time format (default: 'epoch')
            columnnames: Column names to use (auto-detected if None)
            pvformats: Tuple of (position_format, velocity_format) for display

        Returns:
            TimeSeries: AstroPy time series with the data
        """
        # Determine column names based on velocity presence and aux attributes
        if columnnames is None:
            if self.aux:
                # Has aux attributes - use individual position columns
                columnnames = pvhelper._ephemeris_columns_pos_xyz
            elif self.has_velocity:
                # Has velocity - use position and velocity columns
                columnnames = pvhelper._ephemeris_columns
            else:
                # Position only
                columnnames = pvhelper._ephemeris_columns_pos_only
        # Name the columns with the PVT name if available
        if hasattr(self,'name'):
            timenm = columnnames[0]
            columnnames = [self.name + " " + cn for cn in columnnames]
            columnnames[0]=timenm

        # Prepare time array
        if self.time.isscalar:
            tm = astropy.time.Time([self.time.to_string()])
        else:
            tm = self.time

        # Create time series with appropriate data
        if len(columnnames) == 4:  # px, py, pz format
            ts = TimeSeries(time=tm, data=self.position_vector.tovector(), names=columnnames[1:])
        elif len(columnnames) == 2:  # position only
            ts = TimeSeries(time=tm, data=self.position_vector.tovector(), names=columnnames[1:])
        else:  # position and velocity
            ts = TimeSeries(time=tm, data=self.cartesian.tovector(), names=columnnames[1:])

        # Add aux attributes if present
        for key in self.aux:
            ts[key] = self.aux[key].split(' ')

        # Set display formats
        if len(columnnames) == 4:  # px, py, pz format
            ts[columnnames[1]].info.format = pvformats[0]
            ts[columnnames[2]].info.format = pvformats[0]
            ts[columnnames[3]].info.format = pvformats[0]
        elif len(columnnames) == 2:  # position only
            ts[columnnames[1]].info.format = pvformats[0]
        else:  # position and velocity
            ts[columnnames[1]].info.format = pvformats[0]
            ts[columnnames[2]].info.format = pvformats[1]

        # Add elapsed time column
        if elapsed:
            elapsed_times = [dt.quantity_str for dt in np.diff(tm)]
            elapsed_times.insert(0, '')
            ts.add_column(elapsed_times, index=1, name='elapsed')

        # Apply reference time formatting
        if reftime is not None:
            atime.fromtime(ts, reftime=reftime, copy=False)

        return ts

# END class PositionBase()

##################################################
####   PositionT: Position with optional time ####
##################################################

@dataclasses.dataclass
class PositionT(PositionBase):
    """Position in space with optional time (no velocity information).

    Supports lazy conversion between Cartesian and spherical coordinates.
    Can include auxiliary attributes for metadata, flags, etc.
    """

    def __init__(self, time, aux=None, cartesian=None, spherical=None):
        super().__init__(time, aux, cartesian, spherical)

    def _cartesian_to_spherical(self, cartesian):
        """Convert Cartesian state vector to spherical coordinates.

        Args:
            cartesian: Cartesian state vector with position and velocity components

        Returns:
            spherical: Spherical state vector using sph() format
        """
        from astropy.coordinates import CartesianRepresentation, SphericalRepresentation

        pos = cartesian
        if self.isscalar:
            # Create CartesianRepresentation with position
            cart_repr = CartesianRepresentation(x=pos[0], y=pos[1], z=pos[2])
        else:
            # Create CartesianRepresentation with position
            cart_repr = CartesianRepresentation(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2])

        # Convert to spherical representation
        sph_repr = cart_repr.represent_as(SphericalRepresentation)

        # Convert to the format expected by sph() function
        # SphericalRepresentation uses (lon, lat, distance) format
        # which corresponds to (right ascension, declination, distance)
        sphrepr = [sph_repr.lon, sph_repr.lat, sph_repr.distance]
        ret = pvhelper._sphericalpv(sphrepr, None, labels=['rtasc', 'decl', 'distance'])
        return quant.change_units(ret, unit_lookup=units.prefunits)


    def _spherical_to_cartesian(self, spherical):
        """Convert spherical coordinates to Cartesian state vector.

        Args:
            spherical: Spherical state vector in sph() format

        Returns:
            cartesian: Cartesian state vector with position and velocity components
        """
        from astropy.coordinates import CartesianRepresentation, SphericalRepresentation

        # Extract spherical position and velocity from the structured quantity
        rtasc = spherical['rtasc']          # right ascension (longitude)
        decl = spherical['decl']            # declination (latitude)
        distance = spherical['distance']    # radial distance

        # Create SphericalRepresentation with position
        # Note: SphericalRepresentation expects (lon, lat, distance)
        sph_repr = SphericalRepresentation(
            lon=rtasc,
            lat=decl,
            distance=distance
        )

        # Convert to Cartesian representation (this handles both position and velocity)
        cart_repr = sph_repr.represent_as(CartesianRepresentation)
        raise ValueError("This has never been tested")
        # sq = pvhelper._cartesianpv(cart_repr.xyz, unit_lookup=units.prefunits)
        # Create the structured quantity
        return sq

    def copy(self):
        """Create a copy of this PositionT."""
        # Only copy the attribute that's currently defined to avoid unnecessary computation
        if self._cartesian is not None:
            return PositionT(
                time=self.time.copy() if self.time is not None else None,
                aux=self.aux.copy(),
                cartesian=self.cartesian.copy()
            )
        else:
            return PositionT(
                time=self.time.copy() if self.time is not None else None,
                aux=self.aux.copy(),
                spherical=self.spherical.copy()
            )

    def _make_instance(self, time, cartesian, aux):
        """Factory method to create a new PositionT instance."""
        return PositionT(time=time, cartesian=cartesian, aux=aux)

##################################################
####   PositionVelocityT: Full state vector   ####
##################################################

@dataclasses.dataclass
class PositionVelocityT(PositionBase):
    """Orbital state vector as position (Cartesian 3-vector), velocity (Cartesian 3-vector),
    time, and a dictionary of discrete attributes; each field can have multiple rows,
    corresponding to an ephemeris.

    Supports lazy conversion between Cartesian and spherical coordinates.
    """
    def __init__(self, time, aux=None, cartesian=None, spherical=None):
        super().__init__(time, aux, cartesian, spherical)

    @property
    def pv(self):
        """Legacy property name for cartesian coordinates."""
        return self.cartesian

    @pv.setter
    def pv(self, value):
        """Legacy setter for cartesian coordinates."""
        self.cartesian = value

    @property
    def position(self):
        """Extract position-only as a PositionT object."""
        return PositionT(time=self.time, aux=self.aux.copy(), cartesian=self.position_vector)

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
        pos = cartesian[pvhelper._eph_pos]  # position 3-vector
        vel = cartesian[pvhelper._eph_vel]  # velocity 3-vector

        if self.isscalar:
            # Create CartesianRepresentation with position
            cart_repr = CartesianRepresentation(x=pos[0], y=pos[1], z=pos[2])
            # Create CartesianDifferential with velocity
            cart_diff = CartesianDifferential(d_x=vel[0], d_y=vel[1], d_z=vel[2])
        else:
            # Create CartesianRepresentation with position
            cart_repr = CartesianRepresentation(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2])
            # Create CartesianDifferential with velocity
            cart_diff = CartesianDifferential(d_x=vel[:, 0], d_y=vel[:, 1], d_z=vel[:, 2])

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
        ret = pvhelper._sphericalpv(sphrepr, sphrate, labels=['rtasc', 'decl', 'distance'])
        return quant.change_units(ret, unit_lookup=units.prefunits)


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

        raise ValueError("This has never been tested")

        # Extract position and velocity as 3-vectors
        posvel = [cart_repr.x, cart_repr.y, cart_repr.z, \
                  cart_diff.d_x, cart_diff.d_y, cart_diff.d_z]

        # Create the structured quantity
        return pvhelper._cartesianpv(np.array(posvel), None, units.prefunits)

    def copy(self):
        """Create a copy of this PositionVelocityT."""
        # Only copy the attribute that's currently defined to avoid unnecessary computation
        if self._cartesian is not None:
            return PositionVelocityT(
                time=self.time.copy(),
                cartesian=self.cartesian.copy(),
                aux=self.aux.copy()
            )
        else:
            return PositionVelocityT(
                time=self.time.copy(),
                spherical=self.spherical.copy(),
                aux=self.aux.copy()
            )

    def _make_instance(self, time, cartesian, aux):
        """Factory method to create a new PositionVelocityT instance."""
        return PositionVelocityT(time=time, cartesian=cartesian, aux=aux)

    def pvt(self):
        return self

# END class PositionVelocityT

##################################################
####    Make PositionT,  PositionVelocityT    ####
##################################################

def pvtcart(pv, time, specunits=units.prefunits):
    """Define a PositionVelocityT or PositionT by its Cartesian
       components

       pv: A tuple (position, velocity), or list of tuples, or array with 3 (position only) or 6 columns
       time: Any time representation that serve as input to atime.abstime
       specunits: The units to be assigned to the numbers in `pv`

    """
    if isinstance(pv, list):
        cart = quant.vstack([pvhelper._cartesianpv(pv1, True, specunits, units.prefunits) \
                             for pv1 in pv])
    else:
        cart = pvhelper._cartesianpv(pv, None, specunits, units.prefunits)
    if cart.dtype.names and pvhelper._eph_vel in cart.dtype.names:
        return PositionVelocityT(time=atime.abstime(time), cartesian=cart)
    else:
        return PositionT(time=atime.abstime(time), cartesian=cart)

def _pvtattr(object, isscalar):
    """Create PositionVelocityT or PositionT from any object (e.g., ephemeris) that has the time, position, and optionally velocity, properties."""
    hasthing = lambda object, thing: hasattr(object,thing) or (hasattr(object,'colnames') and thing in object.colnames)
    if hasthing(object,pvhelper._eph_time) and hasthing(object,pvhelper._eph_pos):
        if hasthing(object,pvhelper._eph_vel):
            cart = pvhelper._cartesianpv((object[pvhelper._eph_pos], object[pvhelper._eph_vel]), isscalar)
            return PositionVelocityT(time=object[pvhelper._eph_time], cartesian=cart)
        else:
            return PositionT(time=object[pvhelper._eph_time], cartesian=object[pvhelper._eph_pos])
    else:
        raise ValueError('Cannot make a PositionVelocityT or PositionT from this object')

astropy.timeseries.TimeSeries.pvt = lambda self: _pvtattr(self, False)
astropy.table.row.Row.pvt = lambda self: _pvtattr(self, True)

##################################################
#### Legacy aliases for backward compatibility ####
##################################################

# Create aliases for backward compatibility
PVT = PositionVelocityT
