import astro
import numpy as np
import pandas as pd
from astropy.time import Time # "apt" = astropy Time
from datetime import datetime, UTC, timezone
from skyfield import api
skfts = api.load.timescale()
from orekit.pyhelpers import absolutedate_to_datetime, datetime_to_absolutedate
from org.orekit.time import AbsoluteDate # "okad" = Orekit AbsoluteDate

# Convert time in any form to Orekit AbsoluteDate
def to_okad(t):
    if type(t)==Time: # AstroPy
        return datetime_to_absolutedate(t.datetime)
    elif type(t) == np.datetime64: # NumPy
        return datetime_to_absolutedate(pd.Timestamp(t).to_pydatetime())
    elif type(t) == datetime:  # Python
        return datetime_to_absolutedate(t)
    elif type(t) == AbsoluteDate:
        return t

AbsoluteDate.apt = lambda t: Time(absolutedate_to_datetime(t))

# Create the UTC apt at this instant: nowutc()
# Create the Hipparchus AbsoluteDate at this instant: to_ocad(nowutc())
def nowutc(as_okad=False):
    if as_okad:
        return to_okad(nowutc(False))
    else:
        return Time(datetime.now(UTC), scale='utc')

# Convert time to Skyfield time
def to_skftime(t):
    if type(t)==Time: # AstroPy
        return skfts.from_astropy(t)

# To add timezone to datetime
#import pytz
#def utcdt(datetime):
#    return pytz.utc.localize(datetime)
