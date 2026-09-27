import getpass
import sys

import astropy.units as u
from spacetrack import SpaceTrackClient

import tellurion.astro.time as atime

# Example with Sentinel 3A
# isssent = tell.spacetrack_latest([25544, 41335])
# sentst = isssent['SENTINEL 3A']

# NOT AN ACCURATE COMPUTATION OF CARTESIAN POSITION, IT ASSUMES KEPLER ELEMENTS ARE
# OSCULATING:
# sent_badpvt = tell.kepler(sentst.els, sentst.t).pvt() # This no longer works
# A better choice would be to use Orekit to propagate/convert;
# sentgen = tell.prepare(isssent['SENTINEL 3A'], 1*u.day)
# tell.magdiff(sentgen['pvt0'].cartesian, sent_badpvt.cartesian)

class MeanElementSetT:
    """Mean elements, usually defined by fetching from space-track.org."""

    def __init__(self, els, t, tle, model, scdata):
        self.els = els
        self.t = t
        self.tle = tle
        self.model = model
        self.scdata = scdata

    def __repr__(self):
        return (
            f"MeanElementSetT(name={self.scdata.get('name')!r}, "
            f"t={self.t!r}, model={self.model!r})"
        )

    def __eq__(self, other):
        if not isinstance(other, MeanElementSetT):
            return NotImplemented
        return (
            self.tle == other.tle
            and self.model == other.model
            and self.scdata == other.scdata
        )


def satdata(stdict):
    # gp endpoint: EPOCH is a full ISO 8601 string with sub-second precision,
    # e.g. "2026-03-01T12:34:56.789012". No separate EPOCH_MICROSECONDS field.
    epoch = atime.abstime(stdict["EPOCH"])
    orbels = {
        "sma": float(stdict["SEMIMAJOR_AXIS"]) * u.km,
        "ecc": float(stdict["ECCENTRICITY"]) * u.dimensionless_unscaled,
        "inc": float(stdict["INCLINATION"]) * u.deg,
        "raan": float(stdict["RA_OF_ASC_NODE"]) * u.deg,
        "argper": float(stdict["ARG_OF_PERICENTER"]) * u.deg,
        "ma": float(stdict["MEAN_ANOMALY"]) * u.deg,
        "memo": float(stdict["MEAN_MOTION"]) * u.rev / u.day,
        "memod": float(stdict["MEAN_MOTION_DOT"]) * u.rev / (u.day * u.day),
        "memodd": float(stdict["MEAN_MOTION_DDOT"]) * u.rev / (u.day * u.day * u.day),
        "period": float(stdict["PERIOD"]) * u.min,
        "peralt": float(stdict["PERIAPSIS"]) * u.km,  # was PERIGEE in tle_latest
        "apoalt": float(stdict["APOAPSIS"]) * u.km,  # was APOGEE in tle_latest
        "B": 12.7416 * float(stdict["BSTAR"]) * u.m * u.m / u.kg,
    }
    return (orbels, epoch)


def stscdata(stdata):
    return {
        "name": stdata["OBJECT_NAME"],
        "type": stdata["OBJECT_TYPE"],
        "catid": int(stdata["NORAD_CAT_ID"]),  # was OBJECT_NUMBER in tle_latest
        "intldes": stdata["OBJECT_ID"],
    }


def spacetrack_latest(satnums):
    # Get the stclient
    stclient = get_stclient()

    # gp replaces tle_latest; orderby + limit replaces ordinal=1
    stdata = stclient.gp(norad_cat_id=satnums, orderby="epoch desc")
    ret = [
        MeanElementSetT(*satdata(std),
                        [std["TLE_LINE1"], std["TLE_LINE2"]],
                        "SGP4",
                        stscdata(std))
        for std in stdata
    ]
    if type(satnums) is int:
        return ret[0]
    else:
        return dict(zip([el.scdata["name"] for el in ret], ret))


def get_stclient():
    """Prompt for username and password for `space-track.org` and save
    the client globally under `stclient` for use by `spacetrack`
    functions. To avoid the need to enter this information every
    session, define this in your Python startup file:

         import spacetrack
         stclient = spacetrack.SpaceTrackClient(identity="email", password="mypw")
    """
    main_module = sys.modules.get("__main__")
    stclient = getattr(main_module, "stclient", None)
    if stclient:
        return stclient
    else:
        # 1. Prompt the student for credentials securely
        st_identity = input("Enter space-track.org username/email: ")
        st_password = getpass.getpass("Enter space-track.org password: ")
        # 2. Authenticate the client
        stclient = SpaceTrackClient(identity=st_identity, password=st_password)
        del st_password
        # 3. Set the global variable
        setattr(main_module, "stclient", stclient)
        return stclient
