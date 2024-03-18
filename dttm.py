import astro
from astropy.time import Time # "apt" = astropy Time
from datetime import datetime, UTC, timezone
from orekit.pyhelpers import absolutedate_to_datetime, datetime_to_absolutedate
from org.orekit.time import AbsoluteDate # "okad" = Orekit AbsoluteDate

# Convert Astropy Time, absolutedate
Time.okad = lambda t: datetime_to_absolutedate(t.datetime)
AbsoluteDate.apt = lambda t: Time(absolutedate_to_datetime(t))

# Create the UTC apt at this instant: nowutc()
# Create the Hipparchus AbsoluteDate at this instant: nowutc().okad()
def nowutc(as_okad=False):
    if as_okad:
        return nowutc(False).okad()
    else:
        return Time(datetime.now(UTC), scale='utc')
