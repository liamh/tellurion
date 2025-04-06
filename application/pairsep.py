import numpy as np
import numpy.linalg
import scipy.optimize
import astropy.units as u
import tellurion.core as tell
import tellurion.ork as tork

def setmsis(maxdeg, mass, dragarea, dragcoef):
    """Define a forceenv using the NRLMSIS atmospheric drag model with the spacecraft properties defined."""
    return tork.dragforce(tork.setgravity(maxdeg, maxdeg, float(mass.to(u.kg).value)), \
                          'msis', dragcoef=dragcoef, \
                          dragarea= float(dragarea.to(u.m**2).value))

def intrackdeltav(pvt0, mag):
    """Add the in-track velocity and return the new pvt"""
    vel0 = pvt0.pv['velocity']
    uv = vel0/np.linalg.norm(vel0)
    return tell.pvt((pvt0.pv['position'], pvt0.pv['velocity'] + mag*uv, pvt0[1]))

def pairsep(initkep, delay1, dv1, fe1, delay2, dv2, fe2, proptoalt):
    """Propagate the pair of satellites for a time in which the
    unmaneuvered `initkep` will reach an altitude of `proptoalt` after
    one full orbit. The two satellites maneuver with an in-track delta-v after a delay
    from perigee. Returns the separation between the objects.
    """
    initpvt = tork.cartesian(initkep) # Convert Kepler elements to PVT
    period = tork.elementval(initkep, 'period')
    gen = tork.generate(initpvt, 1.5*period, fe1)

    def alttof(tof):
        "Altitude for a given time of flight past the first perigee; set `gen` first"
        pos = tork.propagate(gen, tof[0]*u.s, False).pv['position']
        return np.linalg.norm(pos)-tork.deffe['earthrad']-proptoalt

    # Find the time of flight to the altitude `proptoalt`
    tof = scipy.optimize.fsolve(alttof, (1.1*period).si.value)[0]*u.s

    def maneuver_and_propagate(delay, dv, forceenv):
        pvtatman = intrackdeltav(tork.propagate(gen, delay, False), dv)
        gencan = tork.generate(pvtatman, tof, forceenv)
        return tork.propagate(gencan, tof-delay, False).pv['position']

    pvt1atend = maneuver_and_propagate(delay1, dv1, fe1)
    pvt2atend = maneuver_and_propagate(delay2, dv2, fe2)

    return np.linalg.norm(pvt1atend-pvt2atend)
