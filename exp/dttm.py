"""
Define dates and times and conversions in various packages: Orekit, AstroPy, NumPy, Python
"""

import numpy as np
import pandas as pd
import skyfield.api
import astropy.time
import datetime
import orekit.pyhelpers as pyhelp
import org.orekit.time

def to_okad(t):
    """
    Convert time in any form to Orekit AbsoluteDate (okad)
    """
    if type(t)==astropy.time.Time: # AstroPy
        return pyhelp.datetime_to_absolutedate(t.datetime)
    elif type(t) == np.datetime64: # NumPy
        return pyhelp.datetime_to_absolutedate(pd.Timestamp(t).to_pydatetime())
    elif type(t) == datetime.datetime:  # Python
        return pyhelp.datetime_to_absolutedate(t)
    elif type(t) == org.orekit.time.AbsoluteDate:
        return t

"""
Convert Orekit AbsoluteDate to AstroPy Time
"""
org.orekit.time.AbsoluteDate.apt = lambda t: astropy.time.Time(pyhelp.absolutedate_to_datetime(t))

skfts = skyfield.api.load.timescale()
def to_skftime(t):
    """
    Convert time from AstroPy to Skyfield Time
    """
    if type(t)==astropy.time.Time: # AstroPy
        return skfts.from_astropy(t)

# Create the UTC apt at this instant: nowutc()
# Create the Hipparchus AbsoluteDate at this instant: to_ocad(nowutc())
def nowutc(as_okad=False):
    """
    The time now (in UTC) as Orekit AbsoluteDate (as_okad=True) or AstroPy Time (as_okad=False).
    """
    if as_okad:
        return to_okad(nowutc(False))
    else:
        return astropy.time.Time(datetime.datetime.now(datetime.UTC), scale='utc')

# To add timezone to datetime
#import pytz
#def utcdt(datetime):
#    return pytz.utc.localize(datetime)
