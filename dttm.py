from datetime import datetime, UTC, timezone
from orekit.pyhelpers import absolutedate_to_datetime, datetime_to_absolutedate
from org.orekit.time import AbsoluteDate # "okad" = Orekit AbsoluteDate

# Astropy definitions of time; see https://docs.astropy.org/en/stable/time/index.html
from astropy.time import Time # "apt" = astropy Time
from astropy.units import * # Define time units, e.g. nowutc() + 5*day

# https://docs.astropy.org/en/stable/timeseries/times.html
# from astropy.timeseries import TimeSeries

# Converte Astropy Time, absolutedate
Time.okad = lambda t: datetime_to_absolutedate(t.datetime)
#cannot set 'apt' attribute of immutable type 'datetime.datetime'
#   datetime.apt = lambda t: Time(absolutedate_to_datetime(t))
apt = lambda t: Time(absolutedate_to_datetime(t))

# Create the UTC apt at this instant: nowutc()
# Create the Hipparchus AbsoluteDate at this instant: nowutc().okad()
def nowutc(as_okad=False):
    if as_okad:
        return nowutc(False).okad()
    else:
        return Time(datetime.now(UTC), scale='utc')
