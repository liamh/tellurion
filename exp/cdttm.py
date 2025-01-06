"""
Conversions of datetimes in various packages: Orekit, AstroPy, NumPy, Python
"""

import numpy as np
import pandas as pd
import skyfield.api
import astropy.time
import datetime
import orekit.pyhelpers as pyhelp
import org.orekit.time
import posvel

def okad(t):
    """
    Convert time in any form to Orekit AbsoluteDate (okad)
    """
    if posvel.isdttm(t): # AstroPy
        return pyhelp.datetime_to_absolutedate(t.datetime)
    elif type(t) == np.datetime64: # NumPy
        return pyhelp.datetime_to_absolutedate(pd.Timestamp(t).to_pydatetime())
    elif type(t) == datetime.datetime:  # Python
        return pyhelp.datetime_to_absolutedate(t)
    elif type(t) == org.orekit.time.AbsoluteDate:
        return t

skfts = skyfield.api.load.timescale()
def skftime(t):
    """
    Convert time from dttm to Skyfield Time
    """
    if posvel.isdttm(t):
        return skfts.from_astropy(t)

def dttm(obj):
    '''Convert the object to a dttm as defined by isdttm().'''
    if type(obj) is str:
        return astropy.time.Time(np.datetime64(obj), scale='utc')
    elif type(obj) is org.orekit.time.AbsoluteDate:
        return dttm(pyhelp.absolutedate_to_datetime(obj))
    else:
        return astropy.time.Time(obj)
