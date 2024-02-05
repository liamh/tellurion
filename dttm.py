import datetime
from orekit.pyhelpers import absolutedate_to_datetime
from org.orekit.time import AbsoluteDate, TimeScalesFactory

dtct = {'utc': TimeScalesFactory.getUTC(),
        'hour': 3600.0,
        'day': 86400.0}

def datm(year, month, day, hour=12, minute=0, second=0.0, microsecond=0):
    return(AbsoluteDate(year, month, day, hour, minute, float(second), dtct['utc']))

def nowutc():
    now = datetime.datetime.now(datetime.UTC)
    return(datm(now.year, now.month, now.day, now.hour, now.minute, now.second+1.0e-6*now.microsecond))
