.. _lambert:

************************************
Lambert Transfers (`tell.lambert`)
************************************

Introduction
============

Tellurion provides a Lambert solver through :func:`tellurion.lambert`.
It computes a transfer orbit connecting two position vectors at their
respective epochs.

Basic Usage
===========

Use top-level Tellurion imports and pass two :class:`tellurion.PositionT`
objects with valid epochs::

    import astropy.units as u
    import tellurion as tell

    pvt0 = tell.pvtcart(
        [5740.1326835, 3314.06715, 0.0, -2.7508268, 4.7645718, 5.5016537],
        tell.abstime("2025-01-01T00:00:00"),
    )

    tof = 45.0 * u.minute
    gen = tell.prepare(pvt0, tof, propagator="keplerian")
    pvtf = tell.propagate(gen, [tof], output="pvt")[-1]

    pvt1_lambert, pvt2_lambert = tell.lambert(pvt0.position, pvtf.position)

The function returns two states: the transfer state at the initial epoch and
the transfer state at the final epoch.

Transfer Direction and Revolutions
==================================

Control transfer geometry with:

* ``shortway`` (default ``True``): short-way vs long-way branch.
* ``n_rev`` (default ``0``): number of complete revolutions.

Example::

    pvt_short = tell.lambert(p1, p2, shortway=True, n_rev=0)
    pvt_long = tell.lambert(p1, p2, shortway=False, n_rev=0)

Current Multi-Revolution Limitation
===================================

For ``n_rev > 0``, the underlying Orekit Lambert implementation currently
returns the short-period (low-energy) branch and does not return the
corresponding long-period (high-energy) branch.

Related Tutorial
================

See :doc:`tutorials/lambert1` for a step-by-step notebook workflow.

See Also
========

* :ref:`propagation` for force-model and propagator setup patterns used before
  Lambert targeting workflows.
* :doc:`tutorials/stm1` for sensitivity analysis using propagation Jacobians.

API Reference
=============

.. autofunction:: tellurion.lambert
