.. _posvel:

*******************************************************
Position, Velocity, and Time (`tell.*`)
*******************************************************

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
   :exclude-members: kepler, equinoctial, circular

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
