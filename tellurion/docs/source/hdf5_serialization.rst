.. _hdf5-serialization:

HDF5 Serialization
==================

Tellurion objects can be saved to and loaded from
`HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ files using the
`fsc.hdf5-io <https://fsc-hdf5-io.readthedocs.io/>`_ library, extended
with `astropy-hdf5io <https://github.com/liamh/astropy-hdf5io>`_ for
AstroPy type support.

This works identically for all Tellurion objects — ``PositionVelocityT``,
``ElementSetT``, spacetrack data, or any combination of them.

.. note::
   HDF5 serialization requires the optional ``[hdf5]`` extra. See
   :ref:`hdf5-installation` for installation instructions.

Basic Usage
-----------

Import Tellurion and the ``save``/``load`` functions from ``fsc.hdf5_io``:

.. code-block:: python

   import tellurion as tell
   import astropy.units as u
   from fsc.hdf5_io import save, load

Saving and loading a single object is the same regardless of type:

.. code-block:: python

   # Save
   save(obj, 'filename.hdf5')

   # Load
   obj = load('filename.hdf5')

Examples
--------

Position and Velocity
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import tellurion as tell
   import astropy.units as u
   from fsc.hdf5_io import save, load

   # Create a state vector
   pvt = tell.pvtcart(state_array, time)

   # Save to HDF5
   save(pvt, 'orbit.hdf5')

   # Load back
   loaded_pvt = load('orbit.hdf5')

Orbital Elements
~~~~~~~~~~~~~~~~

.. code-block:: python

   est = tell.kepler(
       tell.allplane({"altper": 400*u.km, "altapo": 400*u.km,
                      "inc": 51.6*u.deg, "argper": 0*u.deg,
                      "raan": 0*u.deg, "ma": 0*u.deg}),
       tell.abstime('2026-01-01 00:00:00'))

   save(est, 'elements.hdf5')
   loaded_est = load('elements.hdf5')

Nested Structures
~~~~~~~~~~~~~~~~~

Dicts and lists of Tellurion and AstroPy objects are saved and
reconstructed correctly:

.. code-block:: python

   from astropy.time import Time

   data = {
       'initial_state': pvt,
       'elements':      est,
       'epoch':         Time('2026-01-01'),
       'altitude':      408 * u.km,
   }

   save(data, 'mission.hdf5')
   loaded = load('mission.hdf5')

   # All types are preserved on load
   print(type(loaded['initial_state']))  # PositionVelocityT
   print(type(loaded['elements']))       # ElementSetT
   print(type(loaded['epoch']))          # astropy.time.Time
   print(type(loaded['altitude']))       # astropy.units.Quantity

Saving to Groups
----------------

For larger projects it is often useful to organise multiple objects
within a single HDF5 file using named groups. The
``save_to_group`` / ``load_from_group`` utilities from
``astropy-hdf5io`` provide this:

.. code-block:: python

   from astropy_hdf5io import save_to_group, load_from_group, print_tree

   # Save different objects to separate groups in one file
   save_to_group(pvt,      'mission.hdf5', 'orbit/initial')
   save_to_group(est,      'mission.hdf5', 'orbit/elements')
   save_to_group(408*u.km, 'mission.hdf5', 'orbit/altitude')

   # Inspect the file structure
   print_tree('mission.hdf5')
   # /
   #   orbit/
   #     initial   [tellurion.PositionVelocityT]
   #     elements  [tellurion.ElementSetT]
   #     altitude  [astropy_hdf5io.Quantity]

   # Load individual items
   loaded_pvt = load_from_group('mission.hdf5', 'orbit/initial')

Saving Nested Structures Recursively
-------------------------------------

For Munch-based or dict-based study structures, ``save_recursive``
walks the hierarchy and saves each leaf object into a matching
group path:

.. code-block:: python

   from munch import Munch
   from astropy_hdf5io import save_recursive, load_recursive

   study = Munch()
   study.init = Munch()
   study.init.pvt  = pvt
   study.init.est  = est
   study.propn = Munch()
   study.propn.pvt = propagated_pvt

   save_recursive(study, 'study.hdf5', 'mission/sat1')

   # Reload as a nested dict (convert back to Munch if desired)
   from munch import munchify
   loaded = munchify(load_recursive('study.hdf5', 'mission/sat1'))

   print(loaded.init.pvt)   # PositionVelocityT
   print(loaded.propn.pvt)  # PositionVelocityT

How It Works
------------

Importing Tellurion automatically registers HDF5 serializers for all
Tellurion types (when the ``[hdf5]`` extra is installed). No additional
imports are required in user code beyond ``from fsc.hdf5_io import save, load``.

Under the hood:

* Tellurion's ``_hdf5`` modules monkey-patch ``to_hdf5()`` methods onto
  each class and register deserializers with ``@subscribe_hdf5``.
* ``astropy-hdf5io`` does the same for all AstroPy types
  (``Quantity``, ``Time``, ``SkyCoord``, ``Table``, etc.).
* ``fsc.hdf5-io`` drives the top-level ``save``/``load``, dispatching
  to these registered serializers automatically.

This means any object that Tellurion or AstroPy knows how to serialise
can be stored in a plain dict or list and saved in one call.