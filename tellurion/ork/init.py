"""
Set up for orbit computations from Orekit
Load with: import orkinit
or from shell
  ipython -m "orkinit" -i
To suppress version information
  OKNOVERPRINT=1 ipython -m "orkinit" -i
"""

import sys, os, astropy #, skyfield
import numpy as np
import orekit
import orekit.pyhelpers as pyhelpers
import pathlib
# from inspect import getmembers # Easy way to see what is defined for an object

# Start Java VM
_vm = orekit.initVM()

if 'OKNOVERPRINT' not in os.environ:  # To suppresss version printing: OKNOVERPRINT=1 ipython
    print ('Python version:',sys.version)
    print ('Java version:',_vm.java_version)
    print ('Numpy version:', np.__version__)
    print ('Astropy version:', astropy.__version__)
    print ('Orekit version:', orekit.VERSION)

# Orekit data path
if 'OREKITDATA' in os.environ:  # set in shell: export OREKITDATA=$(locate orekit-data.zip)
    orkdatapath = os.environ['OREKITDATA']
else:
    thisdir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    envsh = thisdir.parent.absolute().parent.absolute().joinpath("env.sh")
    if pathlib.Path.exists(envsh):
        with open(envsh) as fp:
            lines = [ln for ln in fp]
            orkdatapath = lines[1].split('=')[1].rstrip()
    else:
        orkdatapath = None

# Load the Orekit data file
if orkdatapath == None:
    pyhelpers.setup_orekit_curdir(pathlib.Path("").parent.absolute().parent.absolute()._str)
else:
    pyhelpers.setup_orekit_curdir(orkdatapath)
