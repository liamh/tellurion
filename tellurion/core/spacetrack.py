import itertools

import astropy.units as u

import tellurion.astro.time as atime

# Set stclient using registered space-track.org username and password


# Set stclient using registered space-track.org username and password
# Place these lines with correct username and password
# in ~/.ipython/profile_default/startup/50-spacetrack.py
#   import spacetrack
#   stclient = spacetrack.SpaceTrackClient(identity="myemail@example.com", password="mypw")

# Example with Sentinel 3A
# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sentst = isssent['SENTINEL 3A']
# NOT AN ACCURATE COMPUTATION OF CARTESIAN POSITION, IT ASSUMES KEPLER ELEMENTS ARE
# OSCULATING:
# sent_badpvt = tork.cartesian(tell.kepler(sentst.els, sentst.t))
# A better choice would be to use Orekit to propagate/convert;
# see ork/tle.py for `sent_goodpvt`.
# sentgen = tork.SGP4gen(isssent['SENTINEL 3A'], 1*u.day, {'altitude': 125.0*u.km, 'eclipse': True, 'visibility': []})
# tell.magdiff(sentgen['pvt0'].pv, sent_badpvt.pv)
#

class MeanElementSetT:
    """Mean element set fetched from space-track.org."""

    def __init__(self, els, t, tle, model, scdata):
        self.els = els
        self.t = t
        self.tle = tle
        self.model = model
        self.scdata = scdata

    def __repr__(self):
        return (f"MeanElementSetT(name={self.scdata.get('name')!r}, "
                f"t={self.t!r}, model={self.model!r})")

    def __eq__(self, other):
        if not isinstance(other, MeanElementSetT):
            return NotImplemented
        return (self.tle == other.tle and
                self.model == other.model and
                self.scdata == other.scdata)

def satdata(stdict):
    # gp endpoint: EPOCH is a full ISO 8601 string with sub-second precision,
    # e.g. "2026-03-01T12:34:56.789012". No separate EPOCH_MICROSECONDS field.
    epoch = atime.abstime(stdict["EPOCH"])
    orbels = {"sma": float(stdict["SEMIMAJOR_AXIS"])*u.km,
              "ecc": float(stdict["ECCENTRICITY"])*u.dimensionless_unscaled,
              "inc": float(stdict["INCLINATION"])*u.deg,
              "raan": float(stdict["RA_OF_ASC_NODE"])*u.deg,
              "argper": float(stdict["ARG_OF_PERICENTER"])*u.deg,
              "ma": float(stdict["MEAN_ANOMALY"])*u.deg,
              "memo": float(stdict["MEAN_MOTION"])*u.rev/u.day,
              "memod": float(stdict["MEAN_MOTION_DOT"])*u.rev/(u.day*u.day),
              "memodd": float(stdict["MEAN_MOTION_DDOT"])*u.rev/(u.day*u.day*u.day),
              "period": float(stdict["PERIOD"])*u.min,
              "peralt": float(stdict["PERIAPSIS"])*u.km,   # was PERIGEE in tle_latest
              "apoalt": float(stdict["APOAPSIS"])*u.km,    # was APOGEE in tle_latest
              "B": 12.7416*float(stdict["BSTAR"])*u.m*u.m/u.kg}
    return (orbels, epoch)

def stscdata(stdata):
    return {"name": stdata["OBJECT_NAME"],
            "type": stdata["OBJECT_TYPE"],
            "catid": int(stdata["NORAD_CAT_ID"]),    # was OBJECT_NUMBER in tle_latest
            "intldes": stdata["OBJECT_ID"]}

def spacetrack_latest(stclient, satnums):
    # gp replaces tle_latest; orderby + limit replaces ordinal=1
    stdata = stclient.gp(norad_cat_id=satnums, orderby="epoch desc", limit=1)
    sattle = (stclient.gp(norad_cat_id=satnums, orderby="epoch desc",
                          limit=1, format="tle")).splitlines()
    ret = [MeanElementSetT(*satdata(std), tle, "SGP4", stscdata(std)) \
           for (std, tle) in zip(stdata, itertools.batched(sattle, 2))]
    if type(satnums) is int:
        return ret[0]
    else:
        return dict(zip([el.scdata["name"] for el in ret], ret))
