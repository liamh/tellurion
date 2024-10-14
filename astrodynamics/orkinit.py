# Set up for orbit computations from Orekit
# Load with: from orkinit import *
# or from shell
# ipython -m "orkinit" -i
# OKNOVERPRINT=1 ipython -m "orkinit" -i

# General import
import sys, os, astropy, skyfield
import numpy as np
# Easy way to see what is defined for an object
from inspect import getmembers

# Orekit import and setup
import orekit
import pathlib
vm = orekit.initVM()
if 'OKNOVERPRINT' not in os.environ:  # To suppresss version printing: OKNOVERPRINT=1 ipython
    print ('Python version:',sys.version)
    print ('Java version:',vm.java_version)
    print ('Numpy version:', np.__version__)
    print ('Astropy version:', astropy.__version__)
    print ('Skyfield version:', skyfield.VERSION)
    print ('Orekit version:', orekit.VERSION)

from orekit.pyhelpers import setup_orekit_curdir
# Load the Orekit data file
if 'OREKITDATA' in os.environ:  # set in shell: export OREKITDATA=$(locate orekit-data.zip)
    setup_orekit_curdir(os.environ['OREKITDATA'])
else: # Look in the project directory at the top level
    setup_orekit_curdir(pathlib.Path("").parent.absolute().parent.absolute()._str)

# General math and Orekit utilities
from math import radians, degrees

from org.orekit.utils import Constants
