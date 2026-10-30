import getpass
import sys

import astropy.units as u
import numpy as np
import spacetrack.operators as op
from astropy.time import Time
from spacetrack import SpaceTrackClient

import tellurion.astro.time as atime

# Example with Sentinel 3A
# isssent = tell.spacetrack_latest([25544, 41335])
# sentst = isssent[41335]

# NOT AN ACCURATE COMPUTATION OF CARTESIAN POSITION, IT ASSUMES KEPLER ELEMENTS ARE
# OSCULATING:
# sent_badpvt = tell.kepler(sentst.els, sentst.time).pvt() # This no longer works
# A better choice would be to use Orekit to propagate/convert;
# sentgen = tell.prepare(isssent['SENTINEL 3A'], 1*u.day)
# tell.magdiff(sentgen['pvt0'].cartesian, sent_badpvt.cartesian)

class MeanElementSetT:
    """Mean elements, usually defined by fetching from space-track.org.

    The object is either a single element set or a time series for one
    object.  For a time series, ``time`` is a vector `~astropy.time.Time`,
    each value in ``els`` is a vector `~astropy.units.Quantity`, and
    ``tle`` is ``[line1s, line2s]`` where each is a list of strings.
    For a single set, ``tle`` is ``[line1, line2]``.
    """

    def __init__(self, els, time, tle, model, scdata):
        self.els = els
        self.time = time
        self.tle = tle
        self.model = model
        self.scdata = scdata

    def __repr__(self):
        n = "" if self.time.isscalar else ", n=" + str(self.time.size)
        return (
            f"MeanElementSetT(name={self.scdata.get('name')!r}, "
            f"time={self.time!r}, model={self.model!r}{n})"
        )

    def __eq__(self, other):
        if not isinstance(other, MeanElementSetT):
            return NotImplemented
        return (
            self.tle == other.tle
            and self.model == other.model
            and self.scdata == other.scdata
        )


def _num(stdict, key, unit, scale=1.0):
    val = stdict.get(key)
    return (np.nan if val in (None, "") else scale * float(val)) * unit


def satdata(stdict):
    # gp endpoint: EPOCH is a full ISO 8601 string with sub-second precision,
    # e.g. "2026-03-01T12:34:56.789012". No separate EPOCH_MICROSECONDS field.
    # Missing or null numeric fields become NaN.
    epoch = atime.abstime(stdict["EPOCH"])
    orbels = {
        "sma": _num(stdict, "SEMIMAJOR_AXIS", u.km),
        "ecc": _num(stdict, "ECCENTRICITY", u.dimensionless_unscaled),
        "inc": _num(stdict, "INCLINATION", u.deg),
        "raan": _num(stdict, "RA_OF_ASC_NODE", u.deg),
        "argper": _num(stdict, "ARG_OF_PERICENTER", u.deg),
        "ma": _num(stdict, "MEAN_ANOMALY", u.deg),
        "memo": _num(stdict, "MEAN_MOTION", u.rev / u.day),
        "memod": _num(stdict, "MEAN_MOTION_DOT", u.rev / (u.day * u.day)),
        "memodd": _num(stdict, "MEAN_MOTION_DDOT", u.rev / (u.day * u.day * u.day)),
        "period": _num(stdict, "PERIOD", u.min),
        "peralt": _num(stdict, "PERIAPSIS", u.km),  # was PERIGEE in tle_latest
        "apoalt": _num(stdict, "APOAPSIS", u.km),  # was APOGEE in tle_latest
        "B": _num(stdict, "BSTAR", u.m * u.m / u.kg, 12.7416),
    }
    return (orbels, epoch)


def stscdata(stdata):
    return {
        "name": stdata["OBJECT_NAME"],
        "type": stdata["OBJECT_TYPE"],
        "catid": int(stdata["NORAD_CAT_ID"]),  # was OBJECT_NUMBER in tle_latest
        "intldes": stdata["OBJECT_ID"],
    }


def _meanelset(records, scalar=False):
    """Combine the gp records of one object into a `MeanElementSetT`,
    sorted by increasing epoch.  A single record gives a scalar set if
    ``scalar`` is true, else a length-1 series."""
    records = sorted(records, key=lambda r: r["EPOCH"])
    sats = [satdata(r) for r in records]
    last = records[-1]
    if scalar and len(records) == 1:
        els, time = sats[0]
        tle = [last["TLE_LINE1"], last["TLE_LINE2"]]
    else:
        els = {
            k: u.Quantity([s[0][k].to(sats[0][0][k].unit) for s in sats])
            for k in sats[0][0]
        }
        time = Time([s[1] for s in sats])
        tle = [[r["TLE_LINE1"] for r in records],
               [r["TLE_LINE2"] for r in records]]
    model = last.get("MEAN_ELEMENT_THEORY") or "SGP4"
    return MeanElementSetT(els, time, tle, model, stscdata(last))


def _by_catid(records, scalar=False):
    groups = {}
    for r in records:
        groups.setdefault(int(r["NORAD_CAT_ID"]), []).append(r)
    return {cid: _meanelset(g, scalar) for cid, g in groups.items()}


def spacetrack_query(class_="gp", stclient=None, **predicates):
    """Query space-track.org with any predicates `spacetrack.py` accepts.

    Parameters
    ----------
    class_ : str, optional
        Request class: ``"gp"`` (default, latest element set per
        object), ``"gp_history"`` (all historical element sets),
        or any other class known to `spacetrack.SpaceTrackClient`
        (e.g. ``"satcat"``, ``"decay"``, ``"cdm_public"``).
    stclient : `spacetrack.SpaceTrackClient`, optional
        Client to use; default prompts via `get_stclient`.
    **predicates
        Passed straight to the client: field names such as
        ``norad_cat_id``, ``object_name``, ``epoch``, ``inclination``,
        plus ``orderby``, ``limit``, etc.  Values may use
        `spacetrack.operators`, e.g. ``epoch=op.inclusive_range(a, b)``,
        ``object_name=op.like("STARLINK")``, ``inclination=op.greater_than(98)``.

    Returns
    -------
    dict of {int: `MeanElementSetT`} keyed by NORAD ID for ``gp`` and
    ``gp_history`` (``gp`` gives single sets; ``gp_history`` gives a
    time series per object, ascending in time); otherwise the list of
    raw record dictionaries from space-track.

    Examples
    --------
    >>> import spacetrack.operators as op                    # doctest: +SKIP
    >>> tell.spacetrack_query(object_name=op.like("SENTINEL"),
    ...                       limit=5)                       # doctest: +SKIP
    >>> tell.spacetrack_query("gp_history", norad_cat_id=25544,
    ...     epoch=op.inclusive_range(t0, t1))                # doctest: +SKIP
    """
    if stclient is None:
        stclient = get_stclient()
    records = getattr(stclient, class_)(**predicates)
    if isinstance(records, str):
        raise ValueError(f"space-track class {class_!r} returned non-JSON data")
    if class_ in ("gp", "gp_history"):
        return _by_catid(records, scalar=(class_ == "gp"))
    return list(records)


def spacetrack_latest(satnums, stclient=None):
    """Latest element set for each NORAD ID in ``satnums`` (an int or a
    list of ints).

    Returns
    -------
    dict of {int: `MeanElementSetT`}
        Keyed by NORAD ID, each a single (scalar-time) element set.
    """
    return spacetrack_query("gp", stclient, norad_cat_id=satnums,
                            orderby="epoch desc")


def spacetrack_history(satnums, start, end, stclient=None):
    """Element sets with epochs between ``start`` and ``end`` (inclusive).

    Parameters
    ----------
    satnums : int or list of int
        NORAD catalog IDs.
    start, end : anything accepted by `tellurion.abstime`
        Absolute times bounding the epochs.
    stclient : `spacetrack.SpaceTrackClient`, optional

    Returns
    -------
    dict of {int: `MeanElementSetT`}
        Keyed by NORAD ID, each a time series ascending in epoch
        (``time`` is a vector `~astropy.time.Time`, ``els`` values are
        vector quantities).  IDs with no data in the range are absent.
    """
    rng = op.inclusive_range(*[atime.abstime(t).utc.to_datetime()
                               for t in (start, end)])
    return spacetrack_query("gp_history", stclient, norad_cat_id=satnums,
                            epoch=rng, orderby="epoch asc")


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
