.. _api:

#############
API Reference
#############

Tellurion is intended to be used from a single top-level import::

   import tellurion as tell

The public interface is exposed on ``tell.*``.

Core State and Time
===================

* :func:`tellurion.abstime`
* :func:`tellurion.pvtcart`
* :class:`tellurion.PositionT`
* :class:`tellurion.PositionVelocityT`
  
      
Orbital Elements
================

* :func:`tellurion.kepler`
* :func:`tellurion.equinoctial`
* :func:`tellurion.circular`
* :func:`tellurion.elementval`

Propagation and Force Models
=============================

* :func:`tellurion.prepare`
* :func:`tellurion.propagate`
* :func:`tellurion.setgravity`
* :func:`tellurion.dragforce`

Lambert and Transfer Analysis
=============================

* :func:`tellurion.lambert`

Structured Quantity Utilities
=============================

* :func:`tellurion.make_quantity`
* :func:`tellurion.change_units`
* :func:`tellurion.hstack`
* :func:`tellurion.vstack`
* :func:`tellurion.quantity_to_dict`
* :func:`tellurion.quantity_to_array`
