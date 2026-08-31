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

.. warning::
   Cannot install with pip yet; see below `Installing from Source`_.

The easiest way to install is using pip:

.. code-block:: bash

   pip install tellurion

Installing required packages
----------------------------

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

For the latest development version:

.. code-block:: bash

   git clone https://notavailable-see-liam
   cd tellurion
   pip install -e .

.. note::
   The ``-e`` flag installs in "editable" mode, useful for development.

Development Installation
------------------------

If you want to contribute or modify the code:

.. code-block:: bash

   git clone https://github.com/liamh/your-package.git
   cd your-package
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

Or, if installing from PyPI when available:

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

Orekit/JVM Initialization
-------------------------

Tellurion now initializes Java/Orekit lazily:

* ``import tellurion`` does **not** start the JVM
* the JVM starts automatically on first use of Orekit-backed functionality
* you can also initialize explicitly:

.. code-block:: python

   from tellurion import init_orekit
   init_orekit()

To verify HDF5 support is available:

.. code-block:: python

   from tellurion.core import posvel_hdf5
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
