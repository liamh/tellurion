from datetime import datetime, UTC, timezone
from orekit.pyhelpers import absolutedate_to_datetime, datetime_to_absolutedate
from org.orekit.time import AbsoluteDate

dtct = {'hour': 3600.0,
        'day': 86400.0}

# Python datetime as an ISO8601 string
def isodttm(dttm):
    return(dttm.isoformat())

# Create an AbsoluteDate from the specified datetime components
def datm(year, month, day, hour=12, minute=0, second=0.0):
    return(datetime_to_absolutedate(datetime(year, month, day, hour, minute,
                                             int(second), int(1e6*(second%1.0)),
                                             tzinfo=timezone.utc)))

# The current time as an AbsoluteDate
def nowutc():
    return(datetime_to_absolutedate(datetime.now(UTC)))
