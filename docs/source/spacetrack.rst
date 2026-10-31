.. _spacetrack:

*************************************************************
Space-Track data
*************************************************************

Introduction
============

Tellurion can fetch mean element sets (the data behind two-line
elements) from `space-track.org <https://www.space-track.org>`_ using
the `spacetrack <https://github.com/python-astrodynamics/spacetrack>`_
package (``pip install spacetrack``). You need to `create a free
space-track.org account
<https://www.space-track.org/auth/createAccount>`_ if you do not
already have one. The first call prompts for your username and
password; to avoid that, define ``stclient`` in your Python startup
file::

    import spacetrack
    stclient = spacetrack.SpaceTrackClient(identity="email", password="mypw")

or pass a client explicitly with the ``stclient`` argument of any
function below.

Three functions are provided:

* :func:`~tellurion.spacetrack_latest` -- latest element set by NORAD ID
* :func:`~tellurion.spacetrack_history` -- element sets over a time range
* :func:`~tellurion.spacetrack_query` -- any query that ``spacetrack``
  allows

Results of the first two (and of ``spacetrack_query`` on the ``gp`` and
``gp_history`` classes) are a dictionary keyed by NORAD catalog ID
whose values are :class:`~tellurion.MeanElementSetT`.

Latest element set by NORAD ID
==============================

.. code-block:: python

    import tellurion as tell

    latest = tell.spacetrack_latest([25544, 41335])
    iss = latest[25544]
    iss.scdata["name"]      # 'ISS (ZARYA)'
    iss.time                # epoch, a scalar astropy Time
    iss.els["inc"]          # <Quantity 51.6 deg>
    iss.tle                 # [line1, line2]

By default ``spacetrack_latest`` returns the newest element set however
old it is. Pass ``recent=7`` to accept only epochs within the last 7
days.

Element sets over a time range
==============================

.. code-block:: python

    hist = tell.spacetrack_history(25544, "2026-10-01", "2026-10-07")
    iss = hist[25544]
    len(iss)                # number of element sets
    iss.time                # vector Time, ascending
    iss.els["inc"]          # vector Quantity
    iss.tle                 # [list of line 1s, list of line 2s]

``start`` and ``end`` accept anything :func:`~tellurion.abstime`
accepts. Objects with no data in the range are absent from the result.

To keep only objects still in orbit now, use ``in_orbit=True``. This
first looks up the IDs in ``satcat`` and drops those with a decay date::

    hist = tell.spacetrack_history([25544, 41335], "2026-10-01",
                                   "2026-10-07", in_orbit=True)

General queries
===============

``spacetrack_query(class_="gp", ...)`` passes its keyword arguments to
the ``spacetrack`` client as space-track predicates, so any field of the
request class can be used, together with ``orderby`` and ``limit``.
Comparison, range and pattern values are written with
``spacetrack.operators``.

Find satellites by name
-----------------------

.. code-block:: python

    import spacetrack.operators as op

    sats = tell.spacetrack_query(object_name=op.like("SENTINEL"), limit=10)
    for catid, s in sats.items():
        print(catid, s.scdata["name"])

Select by orbit shape
---------------------

Highly elliptical orbits, perigee altitude below 500 km and apogee
altitude above 35000 km::

    sats = tell.spacetrack_query(
        periapsis=op.less_than(500),
        apoapsis=op.greater_than(35000),
        object_type="PAYLOAD",
        orderby="apoapsis desc",
    )
    for catid, s in sats.items():
        print(catid, s.scdata["name"], s.els["peralt"], s.els["apoalt"])

Select by inclination range
---------------------------

::

    polar = tell.spacetrack_query(inclination=op.inclusive_range(97, 99),
                                  limit=50)

Other request classes
---------------------

For classes other than ``gp`` and ``gp_history`` the result is the list
of raw record dictionaries::

    cat = tell.spacetrack_query("satcat", norad_cat_id=25544, current="Y")
    cat[0]["LAUNCH"]

The ``gp_history`` class can also be queried directly::

    h = tell.spacetrack_query("gp_history", norad_cat_id=25544,
                              epoch=op.inclusive_range("2026-10-01",
                                                       "2026-10-07"))

.. _spacetrack-recent:

Stale element sets and ``recent``
=================================

The ``gp`` class holds every object that has not been marked decayed.
That includes lost and deep-space objects whose last element set is
decades old; for example Giotto (NORAD 15875) has a 1985 epoch.
``spacetrack_query`` therefore restricts ``gp`` queries to epochs
within the last ``recent`` days, **30 by default**::

    tell.spacetrack_query(periapsis=op.less_than(500))             # last 30 days
    tell.spacetrack_query(periapsis=op.less_than(500), recent=7)   # last 7 days
    tell.spacetrack_query(norad_cat_id=15875, recent=None)         # no restriction

``recent`` has no effect for other classes, or if you give an explicit
``epoch`` predicate. Objects in highly elliptical orbits are sometimes
updated less often than every 30 days; widen ``recent`` if results look
sparse.

Time series element sets
========================

A :class:`~tellurion.MeanElementSetT` is either a single element set or
a time series for one object:

.. list-table::
   :header-rows: 1
   :widths: 15, 40, 40

   * - Attribute
     - Single set
     - Time series
   * - ``time``
     - scalar :class:`~astropy.time.Time`
     - vector :class:`~astropy.time.Time`
   * - ``els``
     - dict of scalar quantities
     - dict of vector quantities
   * - ``tle``
     - ``[line1, line2]``
     - ``[list of line1, list of line2]``
   * - ``model``, ``scdata``
     - one value (the newest set's, for a series)
     - same

A series supports ``len()``, indexing and iteration, like a numpy array::

    iss[-1]                 # newest element set, a single set
    iss[:3]                 # first three, a series
    iss[iss.els["inc"] > 51.6 * u.deg]   # boolean mask
    for s in iss:           # one single set per epoch
        print(s.time.iso, s.els["sma"])

A single set is not subscriptable.

Missing or null numeric fields in the space-track record become NaN.

Saving and loading
==================

Element sets of either kind can be saved with the HDF5 support
described in :doc:`hdf5_serialization`.

API Reference
=============

.. autofunction:: tellurion.spacetrack_latest

.. autofunction:: tellurion.spacetrack_history

.. autofunction:: tellurion.spacetrack_query

.. autoclass:: tellurion.MeanElementSetT
   :members:
