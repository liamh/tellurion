.. _posvel:

*******************************************************
Position, Velocity, and Time (`tell.*`)
*******************************************************

.. contents:: Section navigation
   :local:
   :depth: 4

Introduction
============

Tellurion (typically imported as ``import tellurion as tell``) provides classes for representing satellite
orbital states with position, velocity, and time information. It supports:

* Lazy conversion between Cartesian and spherical coordinate systems
* Integration with AstroPy's time and coordinate systems
* Efficient handling of single epochs or time series (ephemerides)
* Auxiliary metadata storage for labels, flags, and other attributes

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

Working with Ephemerides
=========================

Creating Time Series
--------------------

Convert position-velocity states to AstroPy TimeSeries::

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

Auxiliary Attributes
====================

Store metadata with your states::

    pvt = tell.pvtcart(states, None)
    pvt.aux = {
        'satellite_id': 'sat1 sat2 sat3',
        'quality_flag': 'good good suspect'
    }

    # Auxiliary attributes appear in ephemeris
    ts = pvt.ephemeris()
    print(ts['satellite_id'])

Array Conversion
================

Convert to/from numpy arrays::

    # To array (SI units: meters, m/s, MJD)
    arr = pvt.to_array()
    # Shape: (N, 7) for position+velocity+time
    # Or: (N, 4) for position+time only

    # From array
    pvt_new = tell.pvtcart(arr, None)

Example: Complete Workflow
===========================

Here's a complete example using the test data:

.. plot::
   :include-source:

   import numpy as np
   import matplotlib.pyplot as plt
   import tellurion as tell

   # Satellite states (km, km/s, MJD)
   states = np.array([
       [5740132.6835, 3314067.15, 0., -2750.8268, 4764.5718, 5501.6537, 60676.],
       [4581815.8086, 4512263.1753, 1616826.6336, -4891.449, 3141.5916, 5166.4227, 60676.0035],
       [2865499.6272, 5161045.3831, 3036846.5974, -6432.3867, 1140.4153, 4203.822, 60676.0069],
       [801265.5187, 5183536.7281, 4088441.7275, -7187.9582, -990.116, 2736.5135, 60676.0104],
       [-1359960.4135, 4580209.0183, 4646557.5699, -7074.0617, -2989.071, 948.4195, 60676.0139],
       [-3358332.5948, 3427320.1009, 4647312.5703, -6115.1733, -4617.4828, -941.2708, 60676.0174],
       [-4956791.9866, 1865868.1956, 4094285.2497, -4436.2833, -5686.7788, -2706.7532, 60676.0208],
   ])

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

API Reference
=============

Classes
-------

.. autoclass:: tellurion.PositionBase
   :members:

.. autoclass:: tellurion.PositionT
   :members:

.. autoclass:: tellurion.PositionVelocityT
   :members:

.. autoclass:: tellurion.PVT
   :members:

Constructors
------------

.. autofunction:: tellurion.pvtcart

Orekit Conversion Methods on ``PositionVelocityT``
--------------------------------------------------

These instance methods are different from the top-level constructor functions:

* :meth:`tellurion.PositionVelocityT.kepler` converts a Cartesian state to a
  Keplerian element set.
* :meth:`tellurion.PositionVelocityT.equinoctial` converts a Cartesian state to
  an equinoctial element set.
* :meth:`tellurion.PositionVelocityT.circular` converts a Cartesian state to a
  circular element set.
* Top-level :func:`tellurion.kepler`, :func:`tellurion.equinoctial`, and
  :func:`tellurion.circular` construct element sets from element dictionaries.

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
