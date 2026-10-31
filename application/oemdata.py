"""Compatibility import for the Tellurion CCSDS OEM reader."""

import astropy.units as u

from tellurion import oempvt


oemunits = {
    "time": u.second,
    "length": u.km,
    "speed": u.km / u.second,
    "angle": u.degree,
    "angular speed": u.radian / u.second,
    "dimensionless": u.dimensionless_unscaled,
}
