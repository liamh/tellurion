"""
Position, velocity and time sets in AstroPy.

This module provides classes for representing satellite orbital states with
position, velocity, and time information. It supports lazy conversion between
Cartesian and spherical coordinate systems.
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
    """
    Base class for objects with cartesian and/or spherical position information.

    Provides lazy conversion between Cartesian and spherical coordinates,
    with the actual conversion logic implemented by subclasses depending
    on whether velocity information is present.

    Parameters
    ----------
    time : `~astropy.time.Time` or None
        The date and time of the position
    aux : dict, optional
        Discrete attributes (labels, flags, metadata, etc.)
    cartesian : `~astropy.units.Quantity`, optional
        Cartesian coordinates (either position only or position+velocity)
    spherical : `~astropy.units.Quantity`, optional
        Spherical coordinates (right ascension, declination, distance)
    pv : `~astropy.units.Quantity`, optional
        Legacy parameter name for cartesian (for backwards compatibility)

    Attributes
    ----------
    time : `~astropy.time.Time`
        The date and time of the position
    aux : dict
        Discrete attributes (labels, flags, metadata, etc.)
    isscalar : bool
        True if this represents a single epoch rather than multiple

    Notes
    -----
    This is an abstract base class. Use `PositionT` or `PositionVelocityT` instead.

    See Also
    --------
    PositionT : Position without velocity information
    PositionVelocityT : Position with velocity information
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
        """
        Initialize PositionBase with either cartesian or spherical coordinates.

        Parameters
        ----------
        time : `~astropy.time.Time` or None
            Time object or None
        aux : dict, optional
            Dictionary of auxiliary attributes
        cartesian : `~astropy.units.Quantity`, optional
            Cartesian state vector (position + velocity) or array of state vectors
        spherical : `~astropy.units.Quantity`, optional
            Spherical state vector (r, theta, phi, vr, vtheta, vphi)
        pv : `~astropy.units.Quantity`, optional
            Legacy parameter name for cartesian (for backwards compatibility)

        Raises
        ------
        ValueError
            If both cartesian and spherical are specified, or if neither is specified
        """
        self.time = time
        # Handle None time
        self.isscalar = self.time.isscalar if self.time is not None else True
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
        """
        Return the number of epochs in this object.

        Returns
        -------
        int
            Number of epochs

        Raises
        ------
        TypeError
            If called on a scalar object
        """
        if self.isscalar:
            raise TypeError(f"object of type '{type(self).__name__}' has no len()")
        else:
            return len(self.time)

    def __getitem__(self, index):
        """
        Get a single time step or slice.

        Parameters
        ----------
        index : int or slice
            Index or slice to extract

        Returns
        -------
        PositionBase
            A new instance with the selected data

        Raises
        ------
        TypeError
            If called on a scalar object

        Examples
        --------
        >>> pvt = pvtcart(states, None)
        >>> first = pvt[0]  # Get first epoch
        >>> subset = pvt[2:5]  # Get slice
        """
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
        """
        Cartesian coordinates (lazily evaluated from spherical if needed).

        Returns
        -------
        `~astropy.units.Quantity`
            Cartesian coordinates
        """
        if self._cartesian is None:
            self._cartesian = self._spherical_to_cartesian(self._spherical)
        return self._cartesian

    @cartesian.setter
    def cartesian(self, value):
        """
        Set Cartesian coordinates and clear cached spherical.

        Parameters
        ----------
        value : `~astropy.units.Quantity`
            Cartesian coordinates to set
        """
        self._cartesian = value
        self._spherical = None

    @property
    def spherical(self):
        """
        Spherical coordinates (lazily evaluated from Cartesian if needed).

        Returns
        -------
        `~astropy.units.Quantity`
            Spherical coordinates with structured dtype
        """
        if self._spherical is None:
            self._spherical = self._cartesian_to_spherical(self._cartesian)
        return self._spherical

    @spherical.setter
    def spherical(self, value):
        """
        Set spherical coordinates and clear cached Cartesian.

        Parameters
        ----------
        value : `~astropy.units.Quantity`
            Spherical coordinates to set
        """
        self._spherical = value
        self._cartesian = None

    @property
    def position_vector(self):
        """
        Get just the position 3-vector, regardless of whether velocity is present.

        Returns
        -------
        `~astropy.units.Quantity`
            3-vector with position components

        Examples
        --------
        >>> pvt = pvtcart(state, time)
        >>> pos = pvt.position_vector
        >>> pos.shape
        (3,)
        """
        cart = self.cartesian
        # Check if it's a structured quantity with position/velocity fields
        if hasattr(cart, 'dtype') and cart.dtype.names and pvhelper._eph_pos in cart.dtype.names:
            return cart[pvhelper._eph_pos]
        else:
            return cart

    @property
    def has_velocity(self):
        """
        Check if velocity information is present.

        Returns
        -------
        bool
            True if velocity data is available

        Examples
        --------
        >>> pt = PositionT(time=time, cartesian=pos)
        >>> pt.has_velocity
        False
        >>> pvt = PositionVelocityT(time=time, cartesian=state)
        >>> pvt.has_velocity
        True
        """
        cart = self.cartesian
        if hasattr(cart, 'dtype') and hasattr(cart.dtype, 'names') and cart.dtype.names is not None:
            return pvhelper._eph_vel in cart.dtype.names
        return False

    @abc.abstractmethod
    def _cartesian_to_spherical(self, cartesian):
        """
        Convert Cartesian to spherical coordinates.

        Implemented by subclasses to handle position-only vs position+velocity.

        Parameters
        ----------
        cartesian : `~astropy.units.Quantity`
            Cartesian coordinates (format depends on subclass)

        Returns
        -------
        `~astropy.units.Quantity`
            Spherical coordinates in sph() format
        """
        pass

    @abc.abstractmethod
    def _spherical_to_cartesian(self, spherical):
        """
        Convert spherical to Cartesian coordinates.

        Implemented by subclasses to handle position-only vs position+velocity.

        Parameters
        ----------
        spherical : `~astropy.units.Quantity`
            Spherical coordinates in sph() format

        Returns
        -------
        `~astropy.units.Quantity`
            Cartesian coordinates (format depends on subclass)
        """
        pass

    @abc.abstractmethod
    def copy(self):
        """
        Create a copy of this object.

        Implemented by subclasses to return the correct type.

        Returns
        -------
        PositionBase
            A new independent copy
        """
        pass

    @abc.abstractmethod
    def _make_instance(self, time, cartesian, aux):
        """
        Factory method to create a new instance of the correct type.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time object
        cartesian : `~astropy.units.Quantity`
            Cartesian coordinates
        aux : dict
            Auxiliary attributes dictionary

        Returns
        -------
        PositionBase
            New instance of the appropriate subclass
        """
        pass

    def timeorder(self):
        """
        Sort in increasing time order.

        Returns
        -------
        PositionBase
            New instance of same type, sorted by time

        Examples
        --------
        >>> pvt_unordered = pvtcart(states_shuffled, None)
        >>> pvt_ordered = pvt_unordered.timeorder()
        """
        if self.time is None or self.time.isscalar:
            return self.copy()
        sorted_items = sorted(self, key=lambda x: x.time)
        return self._from_sorted_list(sorted_items)

    def _from_sorted_list(self, sorted_list):
        """
        Create instance from a sorted list of single-time instances.

        Parameters
        ----------
        sorted_list : list
            List of single-epoch instances, already sorted

        Returns
        -------
        PositionBase
            Combined instance with all epochs

        Raises
        ------
        ValueError
            If sorted_list is empty
        """
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
        """
        Concatenate rows from another object onto the end of this object.

        Parameters
        ----------
        other : PositionBase or list of PositionBase
        Another instance of the same type, or a list of instances

        Returns
        -------
        PositionBase
        self (modified in place)

        Examples
        --------
        >>> pvt1 = pvtcart(states1, None)
        >>> pvt2 = pvtcart(states2, None)
        >>> pvt1.concatenate(pvt2)
        >>> len(pvt1)
        10
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
            # Use atime.abstime which handles Time concatenation properly
            self.time = atime.abstime([self.time, other.time])

            self.aux = funcy.merge_with(' '.join, self.aux, other.aux) # Merge aux attributes
            self.isscalar = False
        return self

    def merge(self, other):
        """
        Merge with another object and put in time order.

        Parameters
        ----------
        other : PositionBase
            Another instance of the same type

        Returns
        -------
        PositionBase
            New instance of same type, merged and sorted

        Examples
        --------
        >>> pvt_combined = pvt1.merge(pvt2)
        """
        return self.copy().concatenate(other).timeorder()

    def to_array(self, time_format=atime.prefnumabstime):
        """
        Convert to a numpy array using SI units.

        Parameters
        ----------
        time_format : str, optional
            Format for time column (default: 'mjd')

        Returns
        -------
        ndarray
            For `PositionT`: shape (N, 4) with columns [px, py, pz, time]
            For `PositionVelocityT`: shape (N, 7) with columns [px, py, pz, vx, vy, vz, time]
            All physical quantities are in SI units (meters, meters/second).

        Notes
        -----
        Auxiliary attributes are not included in the array output.

        Examples
        --------
        >>> pvt = pvtcart(state, time)
        >>> arr = pvt.to_array()
        >>> arr.shape
        (1, 7)
        >>> arr[0, :3]  # position in meters
        array([5740132683.5, 3314067150. ,          0. ])
        """
        pos = self.position_vector
        has_velocity = self.has_velocity

        if self.time is None:
            # No time information - return as-is without reshaping
            if has_velocity:
                vel = self.cartesian[pvhelper._eph_vel]
                return np.hstack((pos.si.value, vel.si.value))
            else:
                # Just return the position values directly
                return pos.si.value
        else:
            # Has time - ensure 2D output
            if self.time.shape == ():
                time_array = np.array([[self.time.to_value(time_format)]])
                # Reshape position to 2D
                pos_vals = pos.si.value.reshape(1, -1)
            else:
                time_array = self.time.to_value(time_format).reshape(-1, 1)
                pos_vals = pos.si.value.reshape(-1, 3) if pos.si.value.ndim == 1 else pos.si.value

            if has_velocity:
                vel = self.cartesian[pvhelper._eph_vel]
                if self.time.shape == ():
                    vel_vals = vel.si.value.reshape(1, -1)
                else:
                    vel_vals = vel.si.value.reshape(-1, 3) if vel.si.value.ndim == 1 else vel.si.value
                return np.hstack((pos_vals, vel_vals, time_array))
            else:
                return np.hstack((pos_vals, time_array))

    def ephemeris(self, elapsed=True, reftime='epoch', columnnames=None,
                  pvformats=(pvhelper._pos_format, pvhelper._vel_format)):
        """
        Create an AstroPy time series.

        For `PositionT`: Creates a time series with position only.
        For `PositionVelocityT`: Creates a time series with position and velocity.
        Both include any `aux` attributes.

        Parameters
        ----------
        elapsed : bool, optional
            If True (default), add a column with elapsed time from previous step
        reftime : str, optional
            Reference time format (default: 'epoch')
        columnnames : list of str, optional
            Column names to use (auto-detected if None)
        pvformats : tuple of str, optional
            Tuple of (position_format, velocity_format) for display

        Returns
        -------
        `~astropy.timeseries.TimeSeries`
            AstroPy time series with the data

        Examples
        --------
        >>> pvt = pvtcart(SATELLITE_STATES, None)
        >>> ts = pvt.ephemeris()
        >>> print(ts)

        Create ephemeris without elapsed time column:

        >>> ts = pvt.ephemeris(elapsed=False)

        With custom column names:

        >>> ts = pvt.ephemeris(columnnames=['time', 'pos', 'vel'])
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
            tm = astropy.time.Time([self.time.to_value('iso')])
        else:
            tm = self.time

        # Get position data
        pos_data = self.position_vector.tovector()

        # Create time series with appropriate data
        if len(columnnames) == 4:  # px, py, pz, time format
            # Need to split position into individual columns
            ts = TimeSeries(time=tm)
            ts[columnnames[1]] = pos_data[:, 0]
            ts[columnnames[2]] = pos_data[:, 1]
            ts[columnnames[3]] = pos_data[:, 2]
        elif len(columnnames) == 2:  # position only (single column)
            ts = TimeSeries(time=tm, data={columnnames[1]: pos_data})
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
    """
    Position in space with optional time (no velocity information).

    Supports lazy conversion between Cartesian and spherical coordinates.
    Can include auxiliary attributes for metadata, flags, etc.

    Parameters
    ----------
    time : `~astropy.time.Time` or None
        The date and time of the position
    aux : dict, optional
        Dictionary of auxiliary attributes (metadata, labels, flags)
    cartesian : `~astropy.units.Quantity`, optional
        Cartesian position as a 3-vector or N×3 array
    spherical : `~astropy.units.Quantity`, optional
        Spherical position (right ascension, declination, distance)
    pv : `~astropy.units.Quantity`, optional
        Legacy parameter name for cartesian (for backwards compatibility)

    Attributes
    ----------
    cartesian : `~astropy.units.Quantity`
        Cartesian coordinates (lazily computed from spherical if needed)
    spherical : `~astropy.units.Quantity`
        Spherical coordinates (lazily computed from Cartesian if needed)
    position_vector : `~astropy.units.Quantity`
        Just the position 3-vector
    has_velocity : bool
        Always False for PositionT

    Examples
    --------
    Create a position from Cartesian coordinates:

    >>> import astropy.time
    >>> import astropy.units as u
    >>> import numpy as np
    >>> pos_cart = np.array([5740132.6835, 3314067.15, 0.0]) * u.km
    >>> time = astropy.time.Time(60676.0, format='mjd')
    >>> pos = PositionT(time=time, cartesian=pos_cart)

    Access spherical coordinates (computed lazily):

    >>> sph = pos.spherical
    >>> print(sph['rtasc'], sph['decl'], sph['distance'])

    Create from multiple positions:

    >>> positions = np.array([[5740132.6835, 3314067.15, 0.0],
    ...                       [4581815.8086, 4512263.1753, 1616826.6336]]) * u.km
    >>> times = astropy.time.Time([60676.0, 60676.0035], format='mjd')
    >>> pos_series = PositionT(time=times, cartesian=positions)
    >>> len(pos_series)
    2

    See Also
    --------
    PositionVelocityT : Position with velocity information
    pvtcart : Helper function to create position/velocity objects
    """

    def __init__(self, time, aux=None, cartesian=None, spherical=None, pv=None):
        """Initialize PositionT, passing pv parameter to parent."""
        super().__init__(time, aux, cartesian, spherical, pv)

    def _cartesian_to_spherical(self, cartesian):
        """
        Convert Cartesian state vector to spherical coordinates.

        Parameters
        ----------
        cartesian : `~astropy.units.Quantity`
            Cartesian position vector(s)

        Returns
        -------
        `~astropy.units.Quantity`
            Spherical state vector using sph() format with structured dtype
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
        """
        Convert spherical coordinates to Cartesian state vector.

        Parameters
        ----------
        spherical : `~astropy.units.Quantity`
            Spherical state vector in sph() format

        Returns
        -------
        `~astropy.units.Quantity`
            Cartesian position vector(s)
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
        """
        Create a copy of this PositionT.

        Returns
        -------
        PositionT
            A new independent copy

        Examples
        --------
        >>> pos1 = PositionT(time=time, cartesian=pos)
        >>> pos2 = pos1.copy()
        >>> pos2.aux['label'] = 'modified'  # doesn't affect pos1
        """
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
        """
        Factory method to create a new PositionT instance.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time object
        cartesian : `~astropy.units.Quantity`
            Cartesian coordinates
        aux : dict
            Auxiliary attributes dictionary

        Returns
        -------
        PositionT
            New PositionT instance
        """
        return PositionT(time=time, cartesian=cartesian, aux=aux)

##################################################
####   PositionVelocityT: Full state vector   ####
##################################################

@dataclasses.dataclass
class PositionVelocityT(PositionBase):
    """
    Orbital state vector with position, velocity, and time.

    Each field can have multiple rows, corresponding to an ephemeris.
    Supports lazy conversion between Cartesian and spherical coordinates.

    Parameters
    ----------
    time : `~astropy.time.Time`
        The date and time of the state(s)
    aux : dict, optional
        Dictionary of auxiliary attributes (metadata, labels, flags)
    cartesian : `~astropy.units.Quantity`, optional
        Cartesian state vector with 'position' and 'velocity' fields
    spherical : `~astropy.units.Quantity`, optional
        Spherical state vector with position and velocity components
    pv : `~astropy.units.Quantity`, optional
        Legacy parameter name for cartesian (for backwards compatibility)

    Attributes
    ----------
    cartesian : `~astropy.units.Quantity`
        Cartesian state vector (position and velocity)
    spherical : `~astropy.units.Quantity`
        Spherical state vector (lazily computed)
    position_vector : `~astropy.units.Quantity`
        Just the position 3-vector
    position : PositionT
        Extract position-only as a PositionT object
    has_velocity : bool
        Always True for PositionVelocityT
    pv : `~astropy.units.Quantity`
        Legacy alias for cartesian

    Examples
    --------
    Create from a single state vector:

    >>> import numpy as np
    >>> import astropy.time as atime
    >>> state = np.array([5740132.6835, 3314067.15, 0.0,
    ...                   -2750.8268, 4764.5718, 5501.6537])
    >>> time = atime.Time(60676.0, format='mjd')
    >>> pvt = pvtcart(state, time)
    >>> print(pvt.cartesian)

    Create from multiple states with time in array:

    >>> states = np.array([
    ...     [5740132.6835, 3314067.15, 0.0, -2750.8268, 4764.5718, 5501.6537, 60676.0],
    ...     [4581815.8086, 4512263.1753, 1616826.6336, -4891.449, 3141.5916, 5166.4227, 60676.0035]
    ... ])
    >>> pvt_series = pvtcart(states, None)
    >>> ephemeris = pvt_series.ephemeris()

    Extract position only:

    >>> pos = pvt.position
    >>> isinstance(pos, PositionT)
    True

    See Also
    --------
    PositionT : Position without velocity information
    pvtcart : Helper function to create position/velocity objects
    """
    def __init__(self, time, aux=None, cartesian=None, spherical=None, pv=None):
        """Initialize PositionVelocityT, passing pv parameter to parent."""
        super().__init__(time, aux, cartesian, spherical, pv)

    @property
    def pv(self):
        """
        Legacy property name for cartesian coordinates.

        Returns
        -------
        `~astropy.units.Quantity`
            Cartesian state vector (same as cartesian property)
        """
        return self.cartesian

    @pv.setter
    def pv(self, value):
        """
        Legacy setter for cartesian coordinates.

        Parameters
        ----------
        value : `~astropy.units.Quantity`
            Cartesian state vector to set
        """
        self.cartesian = value

    @property
    def position(self):
        """
        Extract position-only as a PositionT object.

        Returns
        -------
        PositionT
            New PositionT object with position but no velocity

        Examples
        --------
        >>> pvt = pvtcart(state, time)
        >>> pos = pvt.position
        >>> pos.has_velocity
        False
        """
        return PositionT(time=self.time, aux=self.aux.copy(), cartesian=self.position_vector)

    def _cartesian_to_spherical(self, cartesian):
        """
        Convert Cartesian state vector to spherical coordinates.

        Parameters
        ----------
        cartesian : `~astropy.units.Quantity`
            Cartesian state vector with position and velocity components

        Returns
        -------
        `~astropy.units.Quantity`
            Spherical state vector using sph() format
        """
        from astropy.coordinates import CartesianRepresentation, SphericalRepresentation
        from astropy.coordinates import CartesianDifferential, SphericalDifferential

        # Extract position and velocity from the structured quantity
        pos = cartesian[pvhelper._eph_pos]  # position 3-vector
        vel = cartesian[pvhelper._eph_vel]  # velocity 3-vector

        # Check actual dimensionality of pos, not just isscalar flag
        if pos.ndim == 1 and pos.shape[0] == 3:
            # 1D array with 3 elements - single position vector
            cart_repr = CartesianRepresentation(x=pos[0], y=pos[1], z=pos[2])
            cart_diff = CartesianDifferential(d_x=vel[0], d_y=vel[1], d_z=vel[2])
        elif pos.ndim == 2:
            # 2D array - multiple positions
            cart_repr = CartesianRepresentation(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2])
            cart_diff = CartesianDifferential(d_x=vel[:, 0], d_y=vel[:, 1], d_z=vel[:, 2])
        else:
            # Fallback - assume scalar
            cart_repr = CartesianRepresentation(x=pos[0], y=pos[1], z=pos[2])
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
        ret = pvhelper._sphericalpv(sphrepr, sphrate, labels=['rtasc', 'decl', 'distance'])
        return quant.change_units(ret, unit_lookup=units.prefunits)


    def _spherical_to_cartesian(self, spherical):
        """
        Convert spherical coordinates to Cartesian state vector.

        Parameters
        ----------
        spherical : `~astropy.units.Quantity`
            Spherical state vector in sph() format

        Returns
        -------
        `~astropy.units.Quantity`
            Cartesian state vector with position and velocity components
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
        """
        Create a copy of this PositionVelocityT.

        Returns
        -------
        PositionVelocityT
            A new independent copy

        Examples
        --------
        >>> pvt1 = pvtcart(state, time)
        >>> pvt2 = pvt1.copy()
        >>> pvt2.aux['label'] = 'modified'  # doesn't affect pvt1
        """
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
        """
        Factory method to create a new PositionVelocityT instance.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time object
        cartesian : `~astropy.units.Quantity`
            Cartesian coordinates
        aux : dict
            Auxiliary attributes dictionary

        Returns
        -------
        PositionVelocityT
            New PositionVelocityT instance
        """
        return PositionVelocityT(time=time, cartesian=cartesian, aux=aux)

    def pvt(self):
        """
        Return self (for API consistency).

        Returns
        -------
        PositionVelocityT
            self

        Notes
        -----
        This method exists for consistency with the `.pvt()` method added to
        AstroPy TimeSeries and Table Row objects.
        """
        return self

# END class PositionVelocityT

##################################################
####    Make PositionT,  PositionVelocityT    ####
##################################################

def pvtcart(pv, time, specunits=units.prefunits):
    """
    Define a PositionVelocityT or PositionT by its Cartesian components.

    Parameters
    ----------
    pv : array-like
        Position and optionally velocity data. Can be:

        - 3-element array: position only [x, y, z] → returns `PositionT`
        - 6-element array: [x, y, z, vx, vy, vz] → returns `PositionVelocityT`
        - N×3 array: multiple positions → returns `PositionT`
        - N×6 array: multiple states → returns `PositionVelocityT`
        - N×4 array: positions with time as last column → returns `PositionT`
        - N×7 array: states with time as last column → returns `PositionVelocityT`
        - list of arrays: each element is processed as above

    time : `~astropy.time.Time` or None
        Time(s) associated with the state(s). If None and `pv` has 4 or 7 columns,
        the last column is interpreted as time in MJD format.
    specunits : dict, optional
        The units to be assigned to the numbers in `pv`. Default is `prefunits`
        which expects km and km/s.

    Returns
    -------
    PositionT or PositionVelocityT
        Appropriate object based on whether velocity is present

    Examples
    --------
    Single state with explicit time:

    >>> state = np.array([5740132.6835, 3314067.15, 0.0,
    ...                   -2750.8268, 4764.5718, 5501.6537])
    >>> time = astropy.time.Time(60676.0, format='mjd')
    >>> pvt = pvtcart(state, time)

    Multiple states with time in array:

    >>> states = np.array([
    ...     [5740132.6835, 3314067.15, 0.0, -2750.8268, 4764.5718, 5501.6537, 60676.0],
    ...     [4581815.8086, 4512263.1753, 1616826.6336, -4891.449, 3141.5916, 5166.4227, 60676.0035]
    ... ])
    >>> pvt = pvtcart(states, None)

    Position only:

    >>> pos = np.array([5740132.6835, 3314067.15, 0.0])
    >>> pt = pvtcart(pos, time)
    >>> pt.has_velocity
    False

    Recreate from propagator output:

    >>> prop = tell.propagate(demoa.propa.gen, prop5m1h, include_init=True, output='pvt')
    >>> new = pvtcart(prop.to_array(), None, tell.siunits)

    See Also
    --------
    PositionT : Position without velocity
    PositionVelocityT : Position with velocity
    """
    if isinstance(pv, list):
        cart = quant.vstack([pvhelper._cartesianpv(pv1, True, specunits, units.prefunits) \
                             for pv1 in pv])
        is_single_state = False
    else:
        is_single_state = False
        if isinstance(pv, np.ndarray):
            width = pv.shape[pv.ndim-1]
            if width==4 or width==7:
                # Handle both 1D and 2D cases
                if pv.ndim == 1:
                    # Single state (1D array with time)
                    is_single_state = True
                    time = atime.from_array(np.array([pv[width-1]]))
                else:
                    # Multiple states with time as last column
                    time = atime.from_array(pv[:,width-1])
        cart = pvhelper._cartesianpv(pv, None, specunits, units.prefunits)

    # Create the appropriate object
    if cart.dtype.names and pvhelper._eph_vel in cart.dtype.names:
        result = PositionVelocityT(time=atime.abstime(time), cartesian=cart)
    else:
        result = PositionT(time=atime.abstime(time), cartesian=cart)

    # Force scalar if input was 1D array with time
    if is_single_state:
        result.isscalar = True

    return result

def _pvtattr(object, isscalar):
    """
    Create PositionVelocityT or PositionT from any object with time and position.

    This is used to add `.pvt()` methods to AstroPy classes.

    Parameters
    ----------
    object : object
        Any object (e.g., ephemeris, TimeSeries, Table Row) that has time,
        position, and optionally velocity properties
    isscalar : bool
        Whether the object represents a single epoch

    Returns
    -------
    PositionVelocityT or PositionT
        Appropriate object based on whether velocity is present

    Raises
    ------
    ValueError
        If the object doesn't have the required time and position attributes

    Notes
    -----
    This function is used internally to monkey-patch `.pvt()` methods onto
    AstroPy's TimeSeries and Table Row classes.
    """
    hasthing = lambda object, thing: hasattr(object,thing) or (hasattr(object,'colnames') and thing in object.colnames)
    if hasthing(object,pvhelper._eph_time) and hasthing(object,pvhelper._eph_pos):
        if hasthing(object,pvhelper._eph_vel):
            cart = pvhelper._cartesianpv((object[pvhelper._eph_pos], object[pvhelper._eph_vel]), isscalar)
            return PositionVelocityT(time=object[pvhelper._eph_time], cartesian=cart)
        else:
            return PositionT(time=object[pvhelper._eph_time], cartesian=object[pvhelper._eph_pos])
    else:
        raise ValueError('Cannot make a PositionVelocityT or PositionT from this object')

# Monkey-patch .pvt() methods onto AstroPy classes
astropy.timeseries.TimeSeries.pvt = lambda self: _pvtattr(self, False)
astropy.table.row.Row.pvt = lambda self: _pvtattr(self, True)

# def vstack(pvts):
#     ret = posvel.PositionVelocityT(time=times, cartesian=\
#                                    quant.vstack(tuple([pvt.cartesian for pvt in pvts])))

##################################################
#### Legacy aliases for backward compatibility ####
##################################################

# Create aliases for backward compatibility
PVT = PositionVelocityT
