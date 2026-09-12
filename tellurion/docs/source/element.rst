.. _element:

*******************************************************
Orbital Elements (`tell.*`)
*******************************************************

Introduction
============

Tellurion (typically imported as ``import tellurion as tell``) provides classes and functions for
representing satellite orbits as sets of scalar orbital elements paired with
an epoch time.  Three classical element types are supported:

* **Keplerian** elements — the traditional six-element set based on conic
  section geometry
* **Equinoctial** elements — a non-singular alternative well-suited to
  near-circular and near-equatorial orbits
* **Circular** elements — a hybrid set that removes the singularity at zero
  eccentricity while retaining inclination and RAAN explicitly

All three types share the same container class :class:`~tellurion.ElementSetT`
and follow the same construction, inspection, and conversion patterns.

Getting Started
===============

Creating a Keplerian Element Set
---------------------------------

The :func:`~tellurion.kepler` factory function accepts a
dictionary of element values.  Any element whose value is not already an
:class:`~astropy.units.Quantity` is given the default unit for its physical
type (see :ref:`units`)::

    import astropy.units as u
    import tellurion as tell

    kep = tell.kepler(
        {"sma": 8000.0*u.km, "ecc": 0.1,
         "inc": 42.0*u.deg, "argper": 66.0*u.deg,
         "raan": 217.4*u.deg, "ma": 7.25*u.deg},
        tell.abstime('2023-09-14T08:30:00'))

    print(kep)               # ElementSetT(elements=..., time=...)
    print(kep.elements)      # structured Quantity
    print(kep.time)          # astropy Time

The epoch argument is optional.  Without it, :func:`~tellurion.kepler`
returns the structured :class:`~astropy.units.Quantity` alone::

    els_only = tell.kepler({"sma": 8000.0*u.km, "ecc": 0.1,
                            "inc": 42.0*u.deg, "argper": 66.0*u.deg,
                            "raan": 217.4*u.deg, "ma": 7.25*u.deg})

Defining Orbits by Altitude
----------------------------

The :func:`~tellurion.allplane` utility converts between the pairs
``(sma, ecc)``, ``(radper, radapo)``, and ``(altper, altapo)`` so that you can
define an orbit in whatever form is most natural::

    # GEO transfer orbit defined by perigee and apogee altitudes
    gto = tell.kepler(
        tell.allplane({"altper": 350*u.km, "altapo": tell.sma(1.0, altitude=True),
                       "inc": 0.0*u.deg, "argper": 120.0*u.deg,
                       "raan": 0.0*u.deg, "ma": 90.0*u.deg}),
        tell.abstime('2026-01-01T05:55:00'))

    # Geosynchronous orbit defined by mean motion
    geo = tell.kepler(
        tell.allplane({"memo": 1.0*u.rev/u.sday, "ecc": 0.0,
                       "inc": 0.0*u.deg, "argper": 120.0*u.deg,
                       "raan": 0.0*u.deg, "ma": 0.0*u.deg}),
        tell.abstime('2026-01-01T20:30:00'))

See :func:`~tellurion.allplane` and :func:`~tellurion.sma`
for the full list of accepted inputs.

Element Set Container
=====================

All element types are stored in an :class:`~tellurion.ElementSetT`
instance — a lightweight class pairing a structured
:class:`~astropy.units.Quantity` with an :class:`~astropy.time.Time` epoch.

Attributes
----------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Attribute
     - Description
   * - ``elements``
     - Structured :class:`~astropy.units.Quantity` whose field names are the
       element names (e.g. ``'sma'``, ``'ecc'``, ``'inc'``).
   * - ``time``
     - :class:`~astropy.time.Time` epoch associated with the element set.

Accessing Individual Elements
------------------------------

Fields of the structured Quantity are accessed by name::

    print(kep.elements['sma'])    # <Quantity 8000. km>
    print(kep.elements['ecc'])    # <Quantity 0.1>
    print(kep.elements['inc'])    # <Quantity 42. deg>

Unpacking
---------

For compatibility with code that expects a two-tuple, :class:`ElementSetT`
supports iteration and indexing::

    els, t = kep          # unpack
    els = kep[0]          # first element → els
    t   = kep[1]          # second element → t

Element Types
=============

Keplerian Elements
------------------

The classical six-element set describing the shape, orientation, and position
of an orbit.

.. list-table::
   :header-rows: 1
   :widths: 15 25 20 40

   * - Name
     - Description
     - Physical type
     - Notes
   * - ``sma``
     - Semi-major axis
     - length
     -
   * - ``ecc``
     - Eccentricity
     - dimensionless
     -
   * - ``inc``
     - Inclination
     - angle
     - Must be in [0°, 180°]
   * - ``argper``
     - Argument of perigee
     - angle
     - ω; normalised to (−180°, 180°]
   * - ``raan``
     - Right ascension of the ascending node
     - angle
     - Ω; normalised to (−180°, 180°]
   * - ``ma``
     - Mean anomaly *(time element, mean variant)*
     - angle
     - M; normalised to (−180°, 180°]
   * - ``ta``
     - True anomaly *(time element, true variant)*
     - angle
     - f; normalised to (−180°, 180°]

Exactly one of ``ma`` or ``ta`` must be present.

Derived quantities available via :func:`~tellurion.elementval`:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Name
     - Description
     - Notes
   * - ``radper``
     - Radius of perigee
     - :math:`a(1-e)`
   * - ``radapo``
     - Radius of apogee
     - :math:`a(1+e)`
   * - ``altper``
     - Altitude of perigee
     - :math:`a(1-e) - R_\oplus`
   * - ``altapo``
     - Altitude of apogee
     - :math:`a(1+e) - R_\oplus`
   * - ``memo``
     - Mean motion
     - :math:`\sqrt{\mu/a^3}`
   * - ``period``
     - Orbital period
     - :math:`2\pi/n`

Equinoctial Elements
--------------------

A non-singular set that avoids the breakdown of Keplerian elements at zero
eccentricity or zero inclination.  The equinoctial elements are related to
the Keplerian elements by:

.. math::

   e_x &= e \cos(\omega + \Omega) \\
   e_y &= e \sin(\omega + \Omega) \\
   h_x &= \tan(i/2) \cos\Omega \\
   h_y &= \tan(i/2) \sin\Omega

.. list-table::
   :header-rows: 1
   :widths: 15 35 20 30

   * - Name
     - Description
     - Physical type
     - Notes
   * - ``sma``
     - Semi-major axis
     - length
     -
   * - ``ex``
     - Equinoctial eccentricity x-component
     - dimensionless
     - :math:`e\cos(\omega+\Omega)`
   * - ``ey``
     - Equinoctial eccentricity y-component
     - dimensionless
     - :math:`e\sin(\omega+\Omega)`
   * - ``hx``
     - Equinoctial inclination x-component
     - dimensionless
     - :math:`\tan(i/2)\cos\Omega`
   * - ``hy``
     - Equinoctial inclination y-component
     - dimensionless
     - :math:`\tan(i/2)\sin\Omega`
   * - ``ml``
     - Mean longitude *(time element, mean variant)*
     - angle
     - :math:`M + \omega + \Omega`
   * - ``tl``
     - True longitude *(time element, true variant)*
     - angle
     - :math:`f + \omega + \Omega`

.. note::
   For a circular orbit ``ex = ey = 0``; for a zero-inclination orbit
   ``hx = hy = 0``.  These are the cases where Keplerian elements are
   singular and equinoctial elements are most valuable.

Circular Elements
-----------------

A hybrid set that removes the eccentricity singularity while keeping
inclination and RAAN in their familiar form.  The eccentricity vector
components are measured relative to the latitude argument:

.. math::

   e_x^c &= e \cos(\omega) \\
   e_y^c &= e \sin(\omega)

.. list-table::
   :header-rows: 1
   :widths: 15 35 20 30

   * - Name
     - Description
     - Physical type
     - Notes
   * - ``sma``
     - Semi-major axis
     - length
     -
   * - ``cex``
     - Circular eccentricity x-component
     - dimensionless
     - :math:`e\cos\omega`
   * - ``cey``
     - Circular eccentricity y-component
     - dimensionless
     - :math:`e\sin\omega`
   * - ``inc``
     - Inclination
     - angle
     -
   * - ``raan``
     - Right ascension of the ascending node
     - angle
     -
   * - ``mla``
     - Mean latitude argument *(time element, mean variant)*
     - angle
     - :math:`M + \omega`
   * - ``tla``
     - True latitude argument *(time element, true variant)*
     - angle
     - :math:`f + \omega`

.. note::
   Circular elements remain singular at zero inclination (the RAAN is
   undefined).  For equatorial near-circular orbits, use equinoctial elements.

Predicates
==========

Three predicate functions identify the type of an element set object.  Each
takes an optional ``est`` keyword that controls whether the argument is
expected to be a full :class:`ElementSetT` (``est=True``, the default) or
just the structured :class:`~astropy.units.Quantity` without an epoch
(``est=False``)::

    tell.iskepels(kep)            # True  — ElementSetT with Keplerian els
    tell.iskepels(kep.els, est=False)  # True  — structured Quantity only

    eq  = pvt.equinoctial()
    tell.isequels(eq)             # True
    tell.iskepels(eq)             # False  — predicates are mutually exclusive

    circ = pvt.circular()
    tell.iscircels(circ)          # True

The three predicates are mutually exclusive: for any :class:`ElementSetT`,
exactly one will return ``True``.

Conversions
===========

Between Element Sets and Cartesian State
-----------------------------------------

The preferred user-facing style is to use instance methods for conversions
between existing Tellurion objects. Use the top-level functions
:func:`~tellurion.kepler`, :func:`~tellurion.equinoctial`,
:func:`~tellurion.circular`, and :func:`~tellurion.pvt` as constructors or
generic adapters::

    # Keplerian → Cartesian → Keplerian
    pvt   = kep.pvt()
    kep2  = pvt.kepler()

    # Cartesian → equinoctial
    eq    = pvt.equinoctial()

    # Cartesian → circular
    circ  = pvt.circular()

    # Equinoctial → Cartesian
    pvt2  = eq.pvt()

    # Circular → Cartesian
    pvt3  = circ.pvt()

.. note::
   The instance conversion methods
   :meth:`~tellurion.ElementSetT.pvt`,
   :meth:`~tellurion.PositionVelocityT.kepler`,
   :meth:`~tellurion.PositionVelocityT.equinoctial`, and
   :meth:`~tellurion.PositionVelocityT.circular` are distinct from the top-level
   constructor functions :func:`~tellurion.kepler`,
   :func:`~tellurion.equinoctial`, and :func:`~tellurion.circular`.
   Use the instance methods when you already have an
   :class:`~tellurion.ElementSetT` or :class:`~tellurion.PositionVelocityT`.
   The instance methods accept an optional ``mean_time_element`` keyword
   (default ``True``) to select between mean and true angle variants.

Reading Individual Element Values
----------------------------------

:func:`~tellurion.elementval` extracts one or more element values
from any orbital state representation and returns them as
:class:`~astropy.units.Quantity` in the preferred units::

    # From a Keplerian ElementSetT
    tell.elementval(kep, 'altper')           # <Quantity 821.86 km>

    # From an equinoctial ElementSetT
    tell.elementval(pvt.equinoctial(), 'sma')

    # Multiple elements at once
    tell.elementval(kep, ['sma', 'ecc', 'inc'])   # returns a list

    # From a Cartesian PVT — dispatches to Keplerian automatically
    tell.elementval(pvt, 'altper')

Time Series of Elements
------------------------

:func:`~tellurion.tselements` builds an AstroPy
:class:`~astropy.timeseries.TimeSeries` of selected elements from a
propagated ephemeris::

    ephem = tell.propagate(gen, times, include_init=True)

    # Keplerian derived quantities
    ts = tell.tselements(ephem, ["altper", "altapo"])

    # Equinoctial components
    ts = tell.tselements(ephem, ["sma", "ex", "ey", "hx", "hy", "ml"])

    # Circular components
    ts = tell.tselements(ephem, ["sma", "cex", "cey", "inc", "raan", "mla"])

The dispatcher inside :func:`~tellurion.tselements` selects the
correct element dictionary automatically based on the names requested.

Example: Complete Workflow
===========================

Here's a complete example converting between all three element types for a
GPS-like semisynchronous orbit::

    import astropy.units as u
    import tellurion as tell

    # Define a semisynchronous orbit (GPS-like)
    kep = tell.kepler(
        tell.allplane({"sma": tell.sma(2.0), "ecc": 0.0,
                       "inc": 55.0*u.deg, "argper": 0.0*u.deg,
                       "raan": 120.0*u.deg, "ma": 77.0*u.deg}),
        tell.abstime('2026-01-01T12:20:00'))

    # Convert to Cartesian and then to the other element types
    pvt  = kep.pvt()
    eq   = pvt.equinoctial()
    circ = pvt.circular()

    # All three give the same semimajor axis
    print(tell.elementval(kep,  'sma').to(u.km))
    print(tell.elementval(eq,   'sma').to(u.km))
    print(tell.elementval(circ, 'sma').to(u.km))

    # Round-trip back to Cartesian
    pvt_from_eq   = eq.pvt()
    pvt_from_circ = circ.pvt()

    # Propagate and extract equinoctial elements over time
    gen  = tell.prepare(pvt, 1*u.day)
    ephem = tell.propagate(gen, np.linspace(5*u.minute, 60*u.minute, 12),
                           include_init=True)
    ts = tell.tselements(ephem, ["sma", "ex", "ey", "ml"])
    print(ts)

Choosing an Element Type
=========================

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Situation
     - Recommended type
   * - General orbit, any eccentricity, any inclination
     - Keplerian
   * - Near-circular orbit (``ecc`` ≲ 0.01)
     - Equinoctial or Circular
   * - Near-equatorial orbit (``inc`` ≲ 1°)
     - Equinoctial
   * - Near-circular *and* near-equatorial
     - Equinoctial
   * - Near-circular but inclined (e.g. ISS, LEO constellations)
     - Circular (retains familiar ``inc`` and ``raan``)
   * - Numerical propagation internal representation
     - Equinoctial (Orekit's default for numerical integrators)

API Reference
=============

Container Class
---------------

.. autoclass:: tellurion.ElementSetT
   :members:

Constructors
------------

.. autofunction:: tellurion.kepler

.. autofunction:: tellurion.equinoctial

.. autofunction:: tellurion.circular

Utilities
---------

.. autofunction:: tellurion.allplane

.. autofunction:: tellurion.sma

.. autofunction:: tellurion.elementval

.. autofunction:: tellurion.tselements

.. seealso::
   :ref:`posvel` — position, velocity, and time representation.

   :ref:`hdf5-serialization` — saving and loading element sets.
