"""
Orbital elements in Orekit
"""

import astropy.units as u
import numpy as np
from astropy.timeseries import TimeSeries

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.orbits import (
    CartesianOrbit,
    CircularOrbit,
    EquinoctialOrbit,
    KeplerianOrbit,
    OrbitType,
    PositionAngleType,
)

from tellurion.astro import units
from tellurion.core import element, posvel
from tellurion.ork import convert, force

###########################################
#### Element values from orbital state ####
###########################################

_eldict = element.sfdict(
    [
        ["sma", "semimajor axis", "length", u.meter, KeplerianOrbit.getA],
        [
            "ecc",
            "eccentricity",
            "dimensionless",
            u.dimensionless_unscaled,
            KeplerianOrbit.getE,
        ],
        ["inc", "inclination", "angle", u.radian, KeplerianOrbit.getI],
        [
            "argper",
            "argument of perigee",
            "angle",
            u.radian,
            KeplerianOrbit.getPerigeeArgument,
        ],
        [
            "radper",
            "radius of perigee",
            "length",
            u.meter,
            lambda kep: kep.getA() * (1.0 - kep.getE()),
        ],
        [
            "radapo",
            "radius of apogee",
            "length",
            u.meter,
            lambda kep: kep.getA() * (1.0 + kep.getE()),
        ],
        [
            "altper",
            "altitude of perigee",
            "length",
            u.meter,
            lambda kep, earthrad: kep.getA() * (1.0 - kep.getE()) - earthrad,
        ],
        [
            "altapo",
            "altitude of apogee",
            "length",
            u.meter,
            lambda kep, earthrad: kep.getA() * (1.0 + kep.getE()) - earthrad,
        ],
        [
            "raan",
            "right ascension of the ascending node",
            "angle",
            u.radian,
            KeplerianOrbit.getRightAscensionOfAscendingNode,
        ],
        ["ta", "true anomaly", "angle", u.radian, KeplerianOrbit.getTrueAnomaly],
        ["ma", "mean anomaly", "angle", u.radian, KeplerianOrbit.getMeanAnomaly],
        [
            "memo",
            "mean motion",
            "angular speed",
            u.radian / u.second,
            KeplerianOrbit.getKeplerianMeanMotion,
        ],
        [
            "period",
            "orbital period",
            "time",
            u.second,
            KeplerianOrbit.getKeplerianPeriod,
        ],
    ]
)

_elphystype = {key: value["phystype"] for key, value in _eldict.items()}


def elementval(
    orbstate, elt, earthrad=force.deffe["earthrad"].si.value, forceenv=force.deffe
):
    """Compute orbital element(s) from any orbital state representation."""
    if element.iskepels(orbstate) or type(orbstate) is KeplerianOrbit:
        ko = _keplerianorbit(orbstate, forceenv=forceenv)
        return element.statefnval(ko, elt, _eldict, earthrad)
    elif element.isequels(orbstate) or type(orbstate) is EquinoctialOrbit:
        eo = _equinoctialorbit(orbstate, forceenv=forceenv)
        return element.statefnval(eo, elt, _eqdict)
    elif element.iscircels(orbstate) or type(orbstate) is CircularOrbit:
        co = _circularorbit(orbstate, forceenv=forceenv)
        return element.statefnval(co, elt, _circdict)
    else:
        # PVT or other Cartesian input: convert to Keplerian
        ko = _keplerianorbit(orbstate, forceenv=forceenv).withCachedPositionAngleType(
            PositionAngleType.MEAN
        )
        return element.statefnval(ko, elt, _eldict, earthrad)


# def elementval(orbstate, elt, earthrad=force.deffe["earthrad"].si.value):
#     """Compute orbital element(s) from any orbital state representation."""
#     if element.iskepels(orbstate) or type(orbstate) is KeplerianOrbit:
#         return element.statefnval(_keplerianorbit(orbstate), elt, _eldict, earthrad)
#     elif element.isequels(orbstate) or type(orbstate) is EquinoctialOrbit:
#         return element.statefnval(_equinoctialorbit(orbstate), elt, _eqdict)
#     elif element.iscircels(orbstate) or type(orbstate) is CircularOrbit:
#         return element.statefnval(_circularorbit(orbstate), elt, _circdict)
#     else:
#         # Fall back to Keplerian conversion (works from PVT etc.)
#         return element.statefnval(_keplerianorbit(orbstate), elt, _eldict, earthrad)


def tselements(ephem, elements, forceenv=force.deffe):
    """Make a time series of selected orbital elements."""
    return TimeSeries(
        time=ephem.time,
        data=[
            dict(zip(elements, elementval(ephrow, elements, forceenv=forceenv)))
            for ephrow in ephem.pvt()
        ],
    )


# def tselements(ephem, elements, forceenv=force.deffe):
#     """Make a time series of selected orbital elements, with angle-type
#     caching for equinoctial orbits to avoid redundant conversions."""

#     # Determine the best angle type to cache based on what's being requested
#     _mean_els = {'ml', 'ma', 'mla'}
#     _true_els = {'tl', 'ta', 'tla'}
#     elset = set(elements)
#     if elset & _mean_els:
#         pat = PositionAngleType.MEAN
#     elif elset & _true_els:
#         pat = PositionAngleType.TRUE
#     else:
#         pat = PositionAngleType.MEAN   # default

#     def _row_elementval(ephrow):
#         orb = _equinoctialorbit(ephrow, forceenv).withCachedPositionAngleType(pat)
#         return element.statefnval(orb, elements, _eqdict)

#     return TimeSeries(time=ephem.time,
#                       data=[dict(zip(elements, _row_elementval(ephrow)))
#                             for ephrow in ephem.pvt()])

# -----------------------------------------
#     Kepler element set
# -----------------------------------------


def _keplerianorbit(
    oes,
    units=(units.prefunits["length"], units.prefunits["angle"]),
    forceenv=force.deffe,
):
    """Make a org.orekit.orbits.KeplerianOrbit from anything"""
    if element.iskepels(oes):
        return _keporb_from_components(oes.elements, oes.time, units, forceenv)
    elif type(oes) is KeplerianOrbit:
        return oes
    elif type(oes) is posvel.PositionVelocityT:
        co = CartesianOrbit(
            convert._tspvc(oes),
            forceenv["celestialframe"],
            forceenv["earthmu"].si.value,
        )
        return OrbitType.KEPLERIAN.convertType(co)
    else:
        raise ValueError("Cannot transform to Keplerian elements")


def _keporb_from_components(
    oes,
    epoch,
    units=(units.prefunits["length"], units.prefunits["angle"]),
    fe=force.deffe,
):
    """Make a org.orekit.orbits.KeplerianOrbit from orbital elements
    as a u.Quantity or Dict"""
    oessi = oes.si.value
    if "ma" in oessi.dtype.names:
        pat = PositionAngleType.MEAN
        anom = float(oessi["ma"])
    elif "ta" in oessi.dtype.names:
        pat = PositionAngleType.TRUE
        anom = float(oessi["ta"])
    else:
        raise ValueError("Time element (ma or ta) required in element set")
    orb = KeplerianOrbit(
        float(oessi["sma"]),
        float(oessi["ecc"]),
        float(oessi["inc"]),
        float(oessi["argper"]),
        float(oessi["raan"]),
        anom,
        pat,
        fe["celestialframe"],
        convert._abstime_to_okad(epoch),
        fe["earthmu"].si.value,
    )
    return orb.withCachedPositionAngleType(pat)


# -----------------------------------------
# Non-Kepler elsets and converters
# -----------------------------------------

# ── Equinoctial element dict ───────────────────────────────────────────────
_eqdict = element.sfdict(
    [
        ["sma", "semimajor axis", "length", u.meter, EquinoctialOrbit.getA],
        [
            "ex",
            "equinoctial ex",
            "dimensionless",
            u.dimensionless_unscaled,
            EquinoctialOrbit.getEquinoctialEx,
        ],
        [
            "ey",
            "equinoctial ey",
            "dimensionless",
            u.dimensionless_unscaled,
            EquinoctialOrbit.getEquinoctialEy,
        ],
        [
            "hx",
            "equinoctial hx",
            "dimensionless",
            u.dimensionless_unscaled,
            EquinoctialOrbit.getHx,
        ],
        [
            "hy",
            "equinoctial hy",
            "dimensionless",
            u.dimensionless_unscaled,
            EquinoctialOrbit.getHy,
        ],
        ["ml", "mean longitude", "angle", u.radian, EquinoctialOrbit.getLM],
        ["tl", "true longitude", "angle", u.radian, EquinoctialOrbit.getLv],
        ["el", "eccentric longitude", "angle", u.radian, EquinoctialOrbit.getLE],
    ]
)

# ── Circular element dict ──────────────────────────────────────────────────
_circdict = element.sfdict(
    [
        ["sma", "semimajor axis", "length", u.meter, CircularOrbit.getA],
        [
            "cex",
            "circular ex",
            "dimensionless",
            u.dimensionless_unscaled,
            CircularOrbit.getCircularEx,
        ],
        [
            "cey",
            "circular ey",
            "dimensionless",
            u.dimensionless_unscaled,
            CircularOrbit.getCircularEy,
        ],
        ["inc", "inclination", "angle", u.radian, CircularOrbit.getI],
        [
            "raan",
            "right ascension of AN",
            "angle",
            u.radian,
            CircularOrbit.getRightAscensionOfAscendingNode,
        ],
        ["mla", "mean latitude argument", "angle", u.radian, CircularOrbit.getAlphaM],
        ["tla", "true latitude argument", "angle", u.radian, CircularOrbit.getAlphaV],
        ["ela", "eccentric lat argument", "angle", u.radian, CircularOrbit.getAlphaE],
    ]
)


def _equinoctialorbit(oes, forceenv=force.deffe):
    """Make an org.orekit.orbits.EquinoctialOrbit from anything."""
    if element.isequels(oes):
        return _eqorb_from_components(oes.elements, oes.time, forceenv)
    elif type(oes) is EquinoctialOrbit:
        return oes
    elif type(oes) is posvel.PositionVelocityT:
        co = CartesianOrbit(
            convert._tspvc(oes),
            forceenv["celestialframe"],
            forceenv["earthmu"].si.value,
        )
        eo = OrbitType.EQUINOCTIAL.convertType(co)
        # Keep co alive to prevent GC while eo is being used
        eo._cartesian_orbit = co
        return eo
    else:
        # Try converting via Keplerian as intermediate
        ko = _keplerianorbit(oes, forceenv=forceenv)
        eo = OrbitType.EQUINOCTIAL.convertType(ko)
        # Keep ko alive
        eo._keplerian_orbit = ko
        return eo


def _eqorb_from_components(oes, epoch, fe=force.deffe):
    """Make an EquinoctialOrbit from a structured element Quantity."""
    oessi = oes.si.value
    if "ml" in oessi.dtype.names:
        pat = PositionAngleType.MEAN
        lon = float(oessi["ml"])
    elif "tl" in oessi.dtype.names:
        pat = PositionAngleType.TRUE
        lon = float(oessi["tl"])
    else:
        raise ValueError("Time element (ml or tl) required in equinoctial set")

    orb = EquinoctialOrbit(
        float(oessi["sma"]),
        float(oessi["ex"]),
        float(oessi["ey"]),
        float(oessi["hx"]),
        float(oessi["hy"]),
        lon,
        pat,
        fe["celestialframe"],
        convert._abstime_to_okad(epoch),
        fe["earthmu"].si.value,
    )

    # Cache the construction angle type so reads of the same type are free
    orb = orb.withCachedPositionAngleType(pat)

    # Defensive check — getCachedPositionAngleType lets us verify this
    if orb.getCachedPositionAngleType() != pat:
        raise RuntimeError(
            f"Orekit cached {orb.getCachedPositionAngleType()} "
            f"but expected {pat}; check Orekit version compatibility"
        )

    return orb


def _circularorbit(oes, forceenv=force.deffe):
    """Make an org.orekit.orbits.CircularOrbit from anything."""
    if element.iscircels(oes):
        return _circorb_from_components(oes.elements, oes.time, forceenv)
    elif type(oes) is CircularOrbit:
        return oes
    elif type(oes) is posvel.PositionVelocityT:
        co = CartesianOrbit(
            convert._tspvc(oes),
            forceenv["celestialframe"],
            forceenv["earthmu"].si.value,
        )
        circ = OrbitType.CIRCULAR.convertType(co)
        circ._cartesian_orbit = co
        return circ
    else:
        ko = _keplerianorbit(oes, forceenv=forceenv)
        circ = OrbitType.CIRCULAR.convertType(ko)
        circ._keplerian_orbit = ko
        return circ


def _circorb_from_components(oes, epoch, fe=force.deffe):
    """Make a CircularOrbit from a structured element Quantity."""
    oessi = oes.si.value
    if "mla" in oessi.dtype.names:
        pat = PositionAngleType.MEAN
        lat = float(oessi["mla"])
    elif "tla" in oessi.dtype.names:
        pat = PositionAngleType.TRUE
        lat = float(oessi["tla"])
    else:
        raise ValueError("Time element (mla or tla) required in circular set")
    orb = CircularOrbit(
        float(oessi["sma"]),
        float(oessi["cex"]),
        float(oessi["cey"]),
        float(oessi["inc"]),
        float(oessi["raan"]),
        lat,
        pat,
        fe["celestialframe"],
        convert._abstime_to_okad(epoch),
        fe["earthmu"].si.value,
    )
    return orb.withCachedPositionAngleType(pat)


###############################
####  Transformations      ####
###############################

# Transformations between Cartesian state vector and Kepler elements


def _kepler(object, forceenv=force.deffe, mean_time_element=True):  # Add prefunits
    """The Kepler element set from the Cartesian PVT or equivalent"""
    co = CartesianOrbit(
        convert._tspvc(object), forceenv["celestialframe"], forceenv["earthmu"].si.value
    )
    ko = OrbitType.KEPLERIAN.convertType(co)
    if mean_time_element:
        elnames = element.kepeltma_names
    else:
        elnames = element.kepeltta_names
    if hasattr(object, "time"):
        dttm = object.time
    kepels = dict(zip(elnames, elementval(ko, elnames)))
    return element.kepler(kepels, dttm)


posvel.PositionVelocityT.kepler = _kepler


def _equinoctial(object, forceenv=force.deffe, mean_time_element=True):
    """The equinoctial element set from the Cartesian PVT or equivalent."""
    co = CartesianOrbit(
        convert._tspvc(object), forceenv["celestialframe"], forceenv["earthmu"].si.value
    )
    eo = OrbitType.EQUINOCTIAL.convertType(co)
    elnames = element.equeltma_names if mean_time_element else element.equeltta_names
    dttm = object.time
    equels = dict(zip(elnames, elementval(eo, elnames)))
    return element.equinoctial(equels, dttm)


def _circular(object, forceenv=force.deffe, mean_time_element=True):
    """The circular element set from the Cartesian PVT or equivalent."""
    co = CartesianOrbit(
        convert._tspvc(object), forceenv["celestialframe"], forceenv["earthmu"].si.value
    )
    circ = OrbitType.CIRCULAR.convertType(co)
    elnames = element.circeltma_names if mean_time_element else element.circeltta_names
    dttm = object.time
    circs = dict(zip(elnames, elementval(circ, elnames)))
    return element.circular(circs, dttm)


posvel.PositionVelocityT.equinoctial = _equinoctial
posvel.PositionVelocityT.circular = _circular


def pvt(object, dttm=None):
    if element.iskepels(object, True):
        return convert._pvt(_keplerianorbit(object))
    elif element.isequels(object, True):
        return convert._pvt(_equinoctialorbit(object))
    elif element.iscircels(object, True):
        return convert._pvt(_circularorbit(object))
    else:
        raise ValueError("Can only transform element sets with an epoch")


def allplane(oesdict, forceenv=force.deffe, unitlookup=units.prefunits):
    """Generate all plane pairs (sma, ecc), (radper, radapo), (altper,
    altapo) from the first or last pairs; additionally, the mean
    motion can be substituted for semimajor axis in the first pair.

    Example 1, convert from altitudes of perigee and apogee to
    semimajor axis and eccentricity
    byalts = tell.kepler({"altper":160*u.km, "altapo":20250*u.km, \
                          "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, \
                          "ma": 0.0*u.deg}, \
                          tell.abstime('2022-02-15T08:30:00'))
    smaecc = tork.allplane(byalts)
    smaecc[0]['sma'] # <Quantity 16583.13646 km>
    smaecc[0]['ecc'] # <Quantity 0.60573583>
    tell.iskepels(byalts) # False
    tell.iskepels(smaecc) # True

    Example 2, convert from semimajor axis and eccentricity to
    altitudes of perigee and apogee
    bysmaecc = {"sma":8000, "ecc":0.1, "inc":45, "argper": 120.0, "raan": 80.0,
               "ma": 0.0}
    alts = tork.allplane(bysmaecc)
    alts[0]['altper'] # <Quantity 821.86354 km>
    alts[0]['altapo'] # <Quantity 2421.86354 km>

    Example 3, define a geosynchronous orbit by mean motion
    geo = tork.allplane({"memo":1.0*u.rev/u.sday, "ecc":0.0*u.dimensionless_unscaled, \
                          "inc":0.0*u.deg, "argper": 120.0*u.deg, "raan": 0.0*u.deg,
                          "ma": 0.0*u.deg})

    Example 4, define a geosynchronous transfer orbit
    gto = tork.allplane({"altper": 350*u.km, "altapo": tork.sma(1.0,True), \
                          "ecc":0.0*u.dimensionless_unscaled, \
                          "inc":0.0*u.deg, "argper": 120.0*u.deg, "raan": 0.0*u.deg,
                          "ma": 0.0*u.deg})

    """
    names = oesdict.keys()
    if ("sma" in names or "memo" in names) and "ecc" in names:  # OR PERIOD IN NAMES
        if "memo" in names:
            smav = sma(oesdict["memo"], False, forceenv, unitlookup)
            new = {
                "sma": smav,
                "radper": smav * (1 - oesdict["ecc"]),
                "radapo": smav * (1 + oesdict["ecc"]),
            }
        else:
            smav = oesdict["sma"]
            memo = np.sqrt(
                forceenv["earthmu"].to(unitlookup["gravconst"]) / smav**3
            ).to(unitlookup["angular speed"], equivalencies=u.dimensionless_angles())
            new = {
                "memo": memo,
                "radper": smav * (1 - oesdict["ecc"]),
                "radapo": smav * (1 + oesdict["ecc"]),
            }
        new["altper"] = new["radper"] - forceenv["earthrad"]
        new["altapo"] = new["radapo"] - forceenv["earthrad"]
    elif "altper" in names and "altapo" in names:
        new = {
            "radper": oesdict["altper"] + forceenv["earthrad"],
            "radapo": oesdict["altapo"] + forceenv["earthrad"],
        }
        new["sma"] = (new["radapo"] + new["radper"]) / 2
        new["ecc"] = (new["radapo"] - new["radper"]) / (new["radapo"] + new["radper"])
    else:
        raise ValueError(
            "Plane must be defined by either (sma, ecc) or (altper, altapo)"
        )
    return oesdict | new


def sma(input, altitude=False, forceenv=force.deffe, unitlookup=units.prefunits):
    """Find the semimajor axis from the mean motion, orbital perioid,
    altitude, or specific energy; if `input` is a number, it is
    assumed to be a mean motion in revolutions/sidereal day.

    Example of geosynchronous satellite semimajor axis
      tork.sma(1.0)
      <Quantity 42164.1696233 km>
    Example of orbital period
      tork.sma(10000*u.s)
      <Quantity 10032.11910363 km>
    Example of specific energy
      tork.sma(-20*(u.km/u.s)**2)

    """
    mu = forceenv["earthmu"]
    if type(input) is u.Quantity:
        quant = input
    else:
        quant = input * 1.0 * u.rev / u.sday
    pdim = u.get_physical_type(quant)
    if pdim == "angular speed":
        s1 = np.cbrt(mu / quant**2)
        s = s1.to(unitlookup["length"], equivalencies=u.dimensionless_angles())
        if altitude:
            return s - forceenv["earthrad"]
        else:
            return s
    elif pdim == "time":
        return sma(u.rev / quant, altitude, forceenv, unitlookup)
    elif pdim == "specific energy":
        return (-mu / (2 * quant)).to(unitlookup["length"])
    elif pdim == "length":
        return quant + forceenv["earthrad"]
    else:
        raise ValueError("Cannot convert quantity to semimajor axis")
