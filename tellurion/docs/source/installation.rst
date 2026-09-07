Installation
============

Requirements
------------

**Required:**

* Orekit (orekit-jpype)
* Python 3.8 or later
* AstroPy 7.0 or later
* NumPy 1.20 or later

**Optional:**

* Matplotlib (for plotting)
* SciPy (for advanced analysis)
* astropy-hdf5io (for HDF5 serialization, see :ref:`hdf5-installation`)

Installing with pip
-------------------

The easiest way to install is using pip:

.. code-block:: bash

   pip install tellurion

Installing dependencies manually (optional)
-------------------------------------------

In most cases, ``pip install tellurion`` is sufficient. The commands below are
only needed if you want to install dependencies individually:

.. code-block:: bash

   pip install jdk4py
   pip install orekit-jpype
   pip install git+https://gitlab.orekit.org/orekit/orekit-data.git
   pip install astropy
   pip install funcy
   pip install pandas
   pip install spacetrack
   pip install geocoder
   pip install munch # for demos
   pip install sphinx sphinx-automodapi sphinx-astropy numpydoc # for documentation

Installing with conda
---------------------

If you use conda/mamba:

.. code-block:: bash

   conda install -c conda-forge tellurion

Installing from Source
----------------------

For the latest development version from source:

.. code-block:: bash

   git clone https://github.com/liamh/tellurion.git
   cd tellurion
   pip install -e .

.. note::
   The ``-e`` flag installs in "editable" mode, useful for development.

Quickstart (Fresh Install)
--------------------------

After a fresh install, run this minimal end-to-end check.

1) Install Tellurion:

.. code-block:: bash

   python -m pip install --upgrade pip
   python -m pip install tellurion

2) Run a quick sanity script:

.. code-block:: python

   import astropy.units as u
   import tellurion as tell

   # Core sanity check (no JVM required)
   t0 = tell.abstime("2025-01-01T00:00:00")
   pvt0 = tell.pvtcart(
       [5740.1326835, 3314.06715, 0.0, -2.7508268, 4.7645718, 5.5016537],
       t0,
   )
   print("Core objects OK:", type(pvt0).__name__)

   # Orekit-backed quickstart (explicitly starts JVM + Orekit data)
   tell.init_orekit()

   gen = tell.prepare(pvt0, 30.0 * u.minute, propagator="keplerian")
   pvtf = tell.propagate(gen, [30.0 * u.minute], output="pvt")[-1]
   print("Propagation OK; final position [m]:", pvtf.position_vector.si.value)

This confirms both the core API and a first propagation run from a clean
environment.

Development Installation
------------------------

If you want to contribute or modify the code:

.. code-block:: bash

   git clone https://github.com/liamh/tellurion.git
   cd tellurion
   pip install -e ".[dev]"

This installs additional development dependencies like pytest and sphinx.

.. _hdf5-installation:

HDF5 Serialization (Optional)
------------------------------

Tellurion can save and load its objects (such as ``PositionVelocityT``,
``ElementSetT``, and space-track data) to HDF5 files. This functionality
requires the optional ``hdf5`` extra, which pulls in
`astropy-hdf5io <https://github.com/liamh/astropy-hdf5io>`_,
`fsc.hdf5-io <https://fsc-hdf5-io.readthedocs.io/>`_, and
`h5py <https://www.h5py.org/>`_.

To install Tellurion with HDF5 support:

.. code-block:: bash

   pip install -e ".[hdf5]"

To install Tellurion with HDF5 support from PyPI:

.. code-block:: bash

   pip install tellurion[hdf5]

Once installed, HDF5 serialization is enabled automatically when you import
Tellurion — no extra import is needed in your code:

.. code-block:: python

   import tellurion as tell
   from fsc.hdf5_io import save, load

   pvt = tell.pvtcart(state, time)
   save(pvt, 'orbit.hdf5')

   loaded = load('orbit.hdf5')

.. note::
   If you try to use HDF5 functions without the ``[hdf5]`` extra installed,
   Tellurion will raise a clear ``ImportError`` with instructions on how to
   install it.

.. seealso::
   :ref:`hdf5-serialization` — saving and loading objects to an HDF5 file.

Verifying Your Installation
----------------------------

To verify the installation worked:

.. code-block:: python

   import tellurion
   print(tellurion.__version__)

Validate with Tests (Source/Dev Install)
----------------------------------------

If you installed from source (especially with ``-e ".[dev]"``), run the
repository tests to validate your environment:

.. code-block:: bash

   pytest tests/

If optional HDF5 dependencies are not installed, run the non-HDF5 subset:

.. code-block:: bash

   pytest tests/ -m "not hdf5"

To run only HDF5 tests:

.. code-block:: bash

   pytest tests/ -m hdf5

Expected outcome guidance:

* A successful installation should complete with no unexpected failures.
* ``skipped`` and ``xfailed`` tests can be normal, depending on optional
  dependencies and platform-specific behavior.
* Exact pass/skip/xfail counts can change over time, so use test status
  (pass/fail) rather than fixed numbers as the primary signal.

Orekit/JVM Initialization
-------------------------

Tellurion now initializes Java/Orekit lazily:

* ``import tellurion`` does **not** initialize the JVM
* the JVM starts automatically on first use of Orekit-backed functionality
  (including instance methods on an existing Cartesian state, e.g.
  ``some_pvt.kepler()``)
* you can also initialize explicitly:

.. code-block:: python

   import tellurion as tell
   tell.init_orekit()

.. note::
   ``tell.init_orekit()`` raises ``RuntimeError`` if Java/Orekit setup fails.
   ``tell.orekit_available()`` is non-raising, but it may initialize Orekit/JVM
   as part of its check.

To verify HDF5 support is available:

.. code-block:: python

   import tellurion as tell
   _ = tell.posvel_hdf5
   print("HDF5 serialization available")

Platform-Specific Notes
-----------------------

Windows
~~~~~~~

.. warning::
   Some dependencies may require Visual Studio Build Tools on Windows.

macOS
~~~~~

If using an M1/M2 Mac, ensure you're using native ARM Python or Rosetta.

Linux
~~~~~

Most Linux distributions work without additional setup.

Troubleshooting
---------------

**ImportError: No module named 'astropy'**

Make sure AstroPy is installed:

.. code-block:: bash

   pip install astropy

**ImportError: HDF5 support requires the optional 'hdf5' extra**

Install the HDF5 optional dependencies:

.. code-block:: bash

   pip install tellurion[hdf5]

**Installation fails with compiler errors**

Try installing pre-built wheels:

.. code-block:: bash

   pip install --only-binary :all: tellurion
