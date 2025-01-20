import os
if os.getenv("OREKITDATA") == None:
    raise ValueError("You need to source env.sh for Orekit to work")

import numpy as np
import astropy.units as u
from astropy.timeseries import TimeSeries
import posvel
import ork.init
import ork.posvel as opv
import ork.element as oel
import ork.force as ofr
import prop
import astro
import cdttm

################ State and elements

p0 = [5740.13268349, 3314.06715   ,    0.]
v0 = [-2.75082684,  4.76457184,  5.50165367]
newyear = cdttm.dttm('2025-01-01T00:00:00')
pvt0 = posvel.pvt((p0,v0,newyear))  # A tuple (Quantity, Time)

# Orekit representations
opvt0 = opv.orkpvt(*pvt0)           # org.orekit.utils.TimeStampedPVCoordinates
ocartorb0 = opvt0.cartesianorbit()  # org.orekit.orbits.CartesianOrbit
okeporb0 = opvt0.keplerianorbit()   # org.orekit.orbits.KeplerianOrbit

# Get the pvt back
pvt_opvt0 = opvt0.pvt()
pvt_ocartorb0 = ocartorb0.pvt()
pvt_okeporb0 = okeporb0.pvt()

################ Propagation two-body

cgen = prop.generate(pvt0, 86400.0) # Use cgen for any propagation up to 1 day

# The example pvt as a CartesianOrbit

proptimes = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
ceph = prop.propagate(cgen, proptimes, True)  # Propagate to each step, and include the initial state in the ephemeris table
ceph_has_pvt0 = posvel.pvt(ceph, 0) == pvt0
pvt15m = posvel.pvt(ceph, 3) # PVT for 15min by index
pvt35m = posvel.pvt(ceph,'2025-01-01 00:35:00') # PVT for 35min by time
pvt45m = posvel.pvt(ceph, 45*u.min) # PVT for 45min by relative time
pvt1h = posvel.pvt(ceph) # PVT at the end of the ephemeris
pvt12h = prop.propagate(cgen, 12*u.hour) # Propagate to a single time, as a PVT

pvtshift = posvel.pvt(pvt12h,1*u.day)

ckep0 = cgen.initialState.orbit.keplerianorbit()

caltperapo = oel.tselements(ceph, ["altper","altapo"]) # TimeTable of altitudes of perigee and apogee
capa = astro.hcat(ceph, caltperapo) # Ephemeris table with additional columns for perige and apogee altitude

# orb1h = opv.orkpvt(*pvt1h).cartesianorbit()
kep1h = oel.kepler(pvt1h)
# pvt1h, opv.pvtork(kep1h) are equal but can't be compared, u.allclose does not work on the pv part
# https://github.com/astropy/astropy/issues/17602, it is a numpy issue https://github.com/numpy/numpy/issues/28104


kep0 = oel.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "raan":217.4, "ma":7.25}, \
                      cdttm.dttm('2023-09-14T08:30:00'))
kgen = prop.generate(kep0, 86400.0)
keph = prop.propagate(kgen, proptimes, True)  # Propagate to each step, and include the initial state in the ephemeris table
# posvel.pvt(keph, 35*u.min)
# fails with KeyError: 'No matches found for key 2023-09-14 09:05:00'
# but posvel.pvt(keph,['2023-09-14 09:05:00']) works

kep20m = oel.kepler(keph[4])

################ Propagation with perturbations

fe4x4 = ofr.setgravity(4,4)
cgen4x4 = prop.generate(pvt0, 86400.0, fe4x4)
ceph4x4 = prop.propagate(cgen4x4, proptimes, True)
ceph4x4_posdiff = posvel.magdiff(ceph4x4['position'], ceph['position'])

# It would be nice to have the ability to assemble table with any columns, appropriately renamed, generalize hcat


############### Earth locations

import astropy.coordinates as coord
import ork.geog as oge
import geog
import geonames

mcd = coord.EarthLocation.of_site('McDonald Observatory')
newyear = cdttm.dttm('2025-01-01T00:00:00')
mcdsv = geog.eciobs(mcd, ceph.time, 'mcdonald sitevec')
cmcd = astro.hcat(ceph, mcdsv) # Ephemeris table with additional column for McDonald site vector
mcdsvork = oge.eciobs(mcd, ceph.time, 'mcdonald sitevec')
cmcdork = astro.hcat(ceph, mcdsvork) # Ephemeris table with additional column for McDonald site vector

# Difference between AstroPy and Orekit
mcdsv_apy = posvel.makepos(geog.eciobs(mcd, newyear))
mcdsv_ork = posvel.makepos(oge.eciobs(mcd, newyear))
mcdsv_apy_ork_dist = np.linalg.norm(mcdsv_ork - mcdsv_apy).si

kickapoo = geog.earthloc(lon='98°45′49.82″W', lat='33°33′08.50″N')
# From geonames.location('carbarn') - integrate earthloc with location()?
carbarn = geog.earthloc(lat=38.87206*u.deg, lon=-77.01748*u.deg, elevation=13.0*u.m)
carbarn_lst_newyear = geog.siderealtime(newyear, carbarn)

############### Earth observations

# A simulated observation from carbarn
ob = geog.azelrange(255*u.deg, 80*u.deg, 1455*u.km, carbarn, newyear)
# ECI Cartesian coordinates
obeci = ob.transform_to(coord.GCRS)
obeci_cart = obeci.cartesian
obeci_radec = obeci.spherical

# An observation from McDonald
obmcd = ob.transform_to(coord.AltAz(location=mcd, obstime=newyear))

# Angles only
# obao = coord.SkyCoord(coord.AltAz(az=63*u.deg, alt=80*u.deg, location=carbarn, obstime=newyear))
# obsloc_carbarn = coord.AltAz(location=carbarn, obstime=newyear)

ptobork = oge.eciobs(carbarn, ob, 'carbarn obs')
ptobapy = geog.eciobs(carbarn, ob, 'carbarn obs')
ob_ork_apy_dist = np.linalg.norm(ptobork.cartesian.xyz - ptobapy.cartesian.xyz).si
