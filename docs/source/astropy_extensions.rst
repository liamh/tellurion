*********************
AstroPy Extensions
*********************

Tellurion extends selected AstroPy classes with convenience methods after::

   import tellurion as tell

TimeSeries Extensions
=====================

These methods are dynamically added to ``astropy.timeseries.TimeSeries``:

* ``TimeSeries.fromtime()``
* ``TimeSeries.components()``

Time Extensions
===============

These methods are dynamically added to ``astropy.time.Time``:

* ``Time.to_array()``

Equivalent top-level function:

.. autofunction:: tellurion.to_array
