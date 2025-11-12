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

Installing with pip
-------------------

.. warning::
   Cannot install with pip yet; see below `Installing from Source`_.

The easiest way to install is using pip:

.. code-block:: bash

   pip install your-package-name

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

   conda install -c conda-forge your-package-name

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

Verifying Your Installation
----------------------------

To verify the installation worked:

.. code-block:: python

   import your_package
   print(your_package.__version__)

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

**Installation fails with compiler errors**

Try installing pre-built wheels:

.. code-block:: bash

   pip install --only-binary :all: your-package-name
