.. _posvel:

*******************************************************
Position, Velocity, and Time
*******************************************************

Introduction
============

Tellurion (typically imported as ``import tellurion as tell``) provides classes for representing satellite
orbital states with position, velocity, and time information. It supports:

* Conversion on demand between Cartesian and spherical coordinate systems
* Efficient handling of single epochs or time series (ephemerides)
* Storage of discrete states and events

Getting Started
===============

Basic Usage
-----------

Creating a position-velocity state from Cartesian coordinates::

    import numpy as np
    import tellurion as tell

    # Single satellite state (position and velocity in km, km/s)
    state = np.array([
        5740.1326835,  3314.06715,  0.0,       # position
        -2.7508268,    4.7645718,   5.5016537  # velocity
    ])
    time = tell.abstime('2025-01-01T00:00:00')

    # Create position-velocity-time object
    pvt = tell.pvtcart(state, time)
    print(pvt.cartesian)
    print(pvt.time)

Working with Multiple Epochs
-----------------------------

Create an ephemeris from multiple states::

    # Array with time as last column (7 columns total)
    states = np.array([
        [5740.1326835,  3314.06715,  0., -2.7508268, 4.7645718, 5.5016537, 60676.0],
        [4581.8158086,  4512.2631753, 1616.8266336, -4.891449, 3.1415916, 5.1664227, 60676.0035],
        [2865.4996272,  5161.0453831, 3036.8465974, -6.4323867, 1.1404153, 4.203822, 60676.0069],
    ])

    # Time is extracted as modified Julian date from the last column
    pvt = tell.pvtcart(states, None)

    # Create AstroPy TimeSeries
    ephemeris = pvt.ephemeris()
    print(ephemeris)

    # Create ephemeris in spherical coordinates
    pvt.ephemeris(coordinate_type='spherical')

Position-Only Objects
---------------------

For position without velocity information::

    import astropy.units as u

    # Just position (3 elements)
    position = np.array([5740132.6835, 3314067.15, 0.0]) * u.km
    time = tell.abstime('2025-01-01T00:00:00')

    pt = tell.pvtcart(position, time)
    print(pt.position_vector)


.. _posvel-frames:

Frame Metadata
--------------

State objects can carry a typed :class:`~tellurion.FrameDescriptor` describing
the coordinate origin and axis convention. This metadata does not resolve an
Orekit frame or perform a coordinate transformation; it makes the state
coordinates explicit for later operations. Existing states default to
``frame=None`` to preserve existing behavior::

    earth_eme2000 = tell.FrameDescriptor(
        origin="Earth",
        axes="EME2000",
        kind=tell.FrameKind.INERTIAL,
    )
    barycentric_icrf = tell.FrameDescriptor(
        origin="Solar System Barycenter",
        axes="ICRF",
        kind=tell.FrameKind.INERTIAL,
    )
    pvt = tell.pvtcart(state, time, frame=earth_eme2000)

Local orbital frames can record their frame kind and, when available, a
reference identifier::

    lvlh_frame = tell.FrameDescriptor(
        origin="reference spacecraft",
        axes="LVLH",
        kind=tell.FrameKind.LOCAL_ORBITAL,
        reference="chief-1",
    )

The reference identifier is descriptive metadata, not the reference
spacecraft trajectory; transformations to or from a local orbital frame
still require that reference state.

Reading CCSDS OEM files
-----------------------

The optional ``oem`` dependency reads CCSDS Orbital Ephemeris Message files
in KVN or XML format. Install it with ``pip install "tellurion[oem]"``; the
dependency is imported only when :func:`tellurion.oempvt` is called. OEM
segment ``CENTER_NAME`` and ``REF_FRAME`` metadata become the PVT's frame
descriptor. OEM state units (km and km/s) are retained::

    pvt = tell.oempvt("trajectory.oem")
    print(pvt.frame)

If consecutive OEM segments change frame, the reader returns a tuple of PVTs,
one per contiguous frame run. This keeps states in unlike frames separate.
Known origins and axes are classified; unknown values are retained in a
custom frame descriptor, so reading does not depend on Orekit supporting that
frame. An unknown descriptor cannot necessarily be transformed by
:func:`tellurion.transform_pvt`.

Querying JPL Horizons
---------------------

The optional Astroquery dependency can fetch geometric state vectors directly
from JPL Horizons. Install it with ``pip install "tellurion[horizons]"``.
The returned PVT uses km and km/s, and its ``frame`` records the requested
center and reference plane. Horizons vector epochs are TDB; ``start_time`` and ``stop_time`` accept
anything :func:`tellurion.abstime` does and are sent as ISO 8601 TDB, and the
returned times display in ISO 8601::

    pvt = tell.horizons(
        target_id="-1024",
        center="@399",
        start_time="2026-04-02T02:00:00",
        stop_time="2026-04-05",
        step="10m",
        refplane="earth",
    )
    print(pvt.frame)

``refplane`` can be ``"earth"``, ``"ecliptic"``, or ``"body"``; its default
is ``"earth"``. ``center="@399"`` identifies Earth, and ``center="@301"``
identifies the Moon. The query uses geometric (uncorrected) vectors, suitable
for representing the target's state rather than an apparent observation.

Orekit Frame Resolution and Transformation
------------------------------------------

Tellurion resolves Earth-centered EME2000 and solar-system-barycentric ICRF
descriptors to Orekit frames. A PVT with a supported descriptor can be
transformed at its epoch, including the frame-origin translation and its
contribution to velocity::

    earth_eme2000 = tell.FrameDescriptor("Earth", "EME2000")
    barycentric_icrf = tell.FrameDescriptor(
        "Solar System Barycenter",
        "ICRF",
    )
    earth_state = tell.pvtcart(state, time, frame=earth_eme2000)
    barycentric_state = tell.transform_pvt(earth_state, barycentric_icrf)
    earth_state_again = tell.transform_pvt(barycentric_state, earth_eme2000)

Use :func:`~tellurion.resolve_orekit_frame` to resolve a descriptor directly.
Unsupported origins, axes, and non-inertial frame kinds are rejected rather
than inferred. Orekit and its data must be initialized before either helper
is used.

Coordinate Systems
==================

The module supports automatic conversion between Cartesian and spherical
coordinates. Conversions are performed lazily - only when accessed.

Cartesian Coordinates
---------------------

Position: :math:`(x, y, z)` in meters

Velocity: :math:`(v_x, v_y, v_z)` in meters/second

Spherical Coordinates
---------------------

Position uses right ascension, declination, and distance:

* **Right ascension** (rtasc): :math:`\alpha`
* **Declination** (decl): :math:`\delta`
* **Distance**: :math:`r`

Velocity uses time derivatives of spherical coordinates.

.. warning::
   Spherical coordinate conversion is currently experimental. The round-trip
   conversion (Cartesian → Spherical → Cartesian) may have numerical precision
   limitations.

.. _ephemeris:

Ephemeris
=========

Creating Time Series
--------------------

Convert position-velocity states to :class:`~astropy.timeseries.TimeSeries`:

    pvt = tell.pvtcart(states, None)
    ts = pvt.ephemeris()

    # Access columns
    print(ts['position'])
    print(ts['velocity'])
    print(ts['elapsed'])  # time between steps

Converting Back from TimeSeries
-------------------------------

AstroPy TimeSeries and Table Row objects have a ``.pvt()`` method added::

    # From TimeSeries
    pvt = ts.pvt()

    # From a single row
    single_state = ts[0].pvt()

.. note::
   The ``.pvt()`` method is automatically added to :class:`~astropy.timeseries.TimeSeries`
   and :class:`~astropy.table.Row` objects when this module is imported.

Manipulating Ephemerides
=========================

Indexing and Slicing
--------------------

Access individual epochs or ranges::

    pvt = tell.pvtcart(states, None)

    # Get first epoch
    first = pvt[0]

    # Get subset
    subset = pvt[2:8]

    # Iterate over epochs
    for epoch in pvt:
        print(epoch.time, epoch.position_vector)

Concatenation and Merging
--------------------------

Combine data from multiple sources::

    pvt1 = tell.pvtcart(states1, None)
    pvt2 = tell.pvtcart(states2, None)

    # Concatenate (in-place)
    pvt1.concatenate(pvt2)

    # Merge and sort by time (returns new object)
    pvt_combined = pvt1.merge(pvt2)

Time Ordering
-------------

Sort ephemerides by time::

    # If data is out of order
    pvt_ordered = pvt.timeorder()

Array Conversion
================

Convert to/from numpy arrays::

    # To array (SI units: meters, m/s, MJD)
    arr = pvt.to_array()
    # Shape: (N, 7) for position+velocity+time
    # Or: (N, 4) for position+time only

    # From array
    pvt_new = tell.pvtcart(arr, None)

Example creating and plotting states
====================================

Here's a complete example using the test data:

.. plot::
   :include-source:

   import numpy as np
   import matplotlib.pyplot as plt
   import tellurion as tell

   # Satellite states (km, km/s, MJD)
   states = np.array([
   [ 5740132.683490001, 3314067.150000001, -0.00000000000000006, -2750.8268400000006, 4764.5718400000005, 5501.65367,60676.],
   [ 4581815.808644368, 4512263.175312268, 1616826.6335879106, -4891.448994499806, 3141.5916389405356, 5166.42266339148, 60676.00347222222],
   [2865499.6272814395, 5161045.383338191, 3036846.5975316903, -6432.386721312811, 1140.4152994517522, 4203.821979628461, 60676.006944444445],
   [801265.5187653155, 5183536.728599629, 4088441.727847057, -7187.958175211956, -990.1159838616379, 2736.5134920974733, 60676.010416666664],
   [-1359960.4135635518, 4580209.016108771, 4646557.568026169, -7074.061713322592, -2989.0709516174024, 948.4194785742409, 60676.01388888889],
   [-3358332.595079403, 3427320.100387147, 4647312.570028321, -6115.173293799234, -4617.482819579691, -941.2707761186798, 60676.01736111111],
   [-4956791.987465892, 1865868.195098941, 4094285.2496156073, -4436.283288001326, -5686.778784978363, -2706.753248705289, 60676.020833333336],
   [-5968503.629846354, 83291.80410506866, 3056384.6323205107, -2243.2158883647585, -6078.38415773859, -4142.427149176573, 60676.024305555555],
   [-6277207.5983626675, -1709253.2344388508, 1658347.07617787, 204.18625334403907, -5753.6241756393, -5084.877825127871, 60676.02777777778],
   [-5848976.973401382, -3301215.6218875158, 65551.89475956949, 2621.9402363866434, -4754.696328239233, -5428.6579241542595, 60676.03125],
   [-4734934.458538764, -4506102.994398553, -1534932.4355012353, 4731.935802527101, -3198.1666168630736, -5135.66143551345, 60676.03472222222],
   [-3065182.5418139943, -5182038.022544307, -2955185.2991340104, 6290.263959665059, -1262.3432189745145, -4238.353274531796, 60676.038194444445],
   [-1034645.9485367368, -5247728.383318046, -4027343.1166760908, 7112.606686790984, 830.5532761846525, -2837.0231062026755, 60676.041666666664]])

   # Create position-velocity-time object
   pvt = tell.pvtcart(states, None)

   # Extract positions
   pos = pvt.position_vector.to('km').value

   # Plot orbit in 3D
   fig = plt.figure(figsize=(8, 8))
   ax = fig.add_subplot(111, projection='3d')
   ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], 'b.-', linewidth=2, markersize=8)
   ax.scatter([0], [0], [0], c='yellow', s=200, marker='o')  # Earth
   ax.set_xlabel('X (km)')
   ax.set_ylabel('Y (km)')
   ax.set_zlabel('Z (km)')
   ax.set_title('Satellite Orbit')
   plt.show()

.. _posvel-relative-frames:

Relative Frames
===============

Tellurion provides relative-state transforms between two spacecraft
states using local orbital frames. There are two common ways of
defining frames; both define the third axis to be the direction of
angular momentum of the orbit. They differ in how the axis in the
orbital plane are defined. In the local-vertical, local-horizontal
(LVLH, or RSW) frame, the first axis aligns with the radial vector of
the reference satellite, and the second is in the general direction of
motion, or "along-track". In the NTW frame, the second axis is
tangential to motion of the reference satellite, in the velocity
vector direction. The terminology "RSW" and "NTW" used comes from
Vallado's "Fundamentals of Astrodynamics and Applications" (Orekit
uses different terms).

The forward functions convert two interial states, or time series of
inertial states, the secondary and primary (or reference) to a (time
series) of relative states of the secondary to the primary.

* :func:`tellurion.lvlh` — LVLH/RSW relative frame
* :func:`tellurion.ntw` — NTW relative frame

The inverse functions convert the relative state(s) and the reference
state(s) to the inertial state(s) of the secondary.

* :func:`tellurion.inv_lvlh` — inverse LVLH/RSW transform (relative -> inertial)
* :func:`tellurion.inv_ntw` — inverse NTW transform (relative -> inertial)

These functions work with single states or aligned ephemerides; that
is, the times must match between secondary and reference states, and
have the same length. The functions have common arguments.

* ``relsc`` is the secondary state (or series).
* ``refsc`` is the reference or primary inertial state (or series) that
  defines frame origin and orientation.
* ``unitlookup`` controls output units for forward transforms and defaults to
  preferred units (pass ``tell.siunits`` for SI).

LVLH (RSW) frame
----------------

The LVLH frame is centered on the reference spacecraft with three axes,

* **R** — radial, from central body toward the reference spacecraft,
* **S** — along-track, completing the right-handed set,
* **W** — cross-track, along orbital angular momentum.

Use :func:`tellurion.lvlh` for inertial -> LVLH relative coordinates and
:func:`tellurion.inv_lvlh` for LVLH -> inertial reconstruction.

.. code-block:: python

    import tellurion as tell

    t0 = tell.abstime('2025-01-01T00:00:00')
    primary = tell.pvtcart([5740.1326835, 3314.06715, 0.0,
                          -2.7508268, 4.7645718, 5.5016537], t0)
    secondary = tell.pvtcart([5739.0, 3315.2, 0.4,
                           -2.75, 4.76, 5.50], t0)

    rel_lvlh = tell.lvlh(secondary, primary)              # preferred units
    rel_lvlh_si = tell.lvlh(secondary, primary, tell.siunits)

    secondary_reconstructed = tell.inv_lvlh(rel_lvlh, primary)

NTW frame
---------

The NTW frame is centered on the reference spacecraft with axes:

* **N** — in the orbital plane, normal to the velocity vector
* **T** — tangential, along the velocity vector (in-track)
* **W** — cross-track, along orbital angular momentum

Use :func:`tellurion.ntw` for inertial -> NTW relative coordinates and
:func:`tellurion.inv_ntw` for NTW -> inertial reconstruction.

.. code-block:: python

    import tellurion as tell

    t0 = tell.abstime('2025-01-01T00:00:00')
    primary = tell.pvtcart([5740.1326835, 3314.06715, 0.0,
                          -2.7508268, 4.7645718, 5.5016537], t0)
    secondary = tell.pvtcart([5739.0, 3315.2, 0.4,
                           -2.75, 4.76, 5.50], t0)

    rel_ntw = tell.ntw(secondary, primary, tell.siunits)

    secondary_reconstructed = tell.inv_ntw(rel_ntw, primary)


API Reference
=============

Classes
-------

.. autoclass:: tellurion.FrameDescriptor
   :members:

.. autoclass:: tellurion.FrameKind
   :members:

.. autoclass:: tellurion.PositionBase
   :members:

.. autoclass:: tellurion.PositionT
   :members:

.. autoclass:: tellurion.PositionVelocityT
   :members:
   :exclude-members: kepler, equinoctial, circular

.. autoclass:: tellurion.PVT
   :members:      

Constructors
------------

.. autofunction:: tellurion.pvtcart

.. autofunction:: tellurion.oempvt

.. autofunction:: tellurion.horizons

.. autofunction:: tellurion.resolve_orekit_frame

.. autofunction:: tellurion.transform_pvt

Helpers
-------

.. autofunction:: tellurion.abstime

Conversion to elements
----------------------

These instance methods convert the orbital state to one of the element set representations :class:`ElementSetT <tellurion.ElementSetT>`.

* :meth:`tellurion.PositionVelocityT.kepler` converts an orbital state to a
  Keplerian element set.
* :meth:`tellurion.PositionVelocityT.equinoctial` converts an orbital state to
  an equinoctial element set.
* :meth:`tellurion.PositionVelocityT.circular` converts an orbital state to a
  circular element set.

Note that the functions :func:`tellurion.kepler`, :func:`tellurion.equinoctial`, and :func:`tellurion.circular` are used to construct element sets from element dictionaries.

``PositionVelocityT.kepler``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tellurion.PositionVelocityT.kepler

``PositionVelocityT.equinoctial``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tellurion.PositionVelocityT.equinoctial

``PositionVelocityT.circular``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tellurion.PositionVelocityT.circular

.. seealso::
   :ref:`hdf5-serialization` — saving and loading ``PositionVelocityT`` objects.
