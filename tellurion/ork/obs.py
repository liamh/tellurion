import numpy as np
from org.orekit.frames import TopocentricFrame
from tellurion.astro import units
from tellurion.core import posvel
from tellurion.core import obs
from tellurion.ork import force
from tellurion.ork import convert
from tellurion.ork import geog

###############################
#### Observation from PVT  ####
###############################

# Make use of element.statefnval
# element.statefnval(pvt, ['azim', 'elev'], obsdict, location)

def _topoframe(loc, name="topo loc", forceenv=force.deffe):
    return TopocentricFrame(forceenv['earth'], geog.geodpt(loc), name)

# Example demoa.prop.eclipse.ephempvt
# ob1 = tork.aer(demoa.init.pvt, eloc.mcd)
def aer(pvt, location, name="topo loc", forceenv=force.deffe):
    '''The az-el-range for a PVT'''
    tspvc = convert._tspvc(pvt)
    oktime = tspvc.getDate()
    okpos = tspvc.getPosition()
    tf = _topoframe(location, name, forceenv)
    cf = forceenv['celestialframe']
    azm = tf.getAzimuth(okpos, cf, oktime)*units.orkunits["angle"]
    elv = tf.getElevation(okpos, cf, oktime)*units.orkunits["angle"]
    rng = tf.getRange(okpos, cf, oktime)*units.orkunits["length"]
    return obs.azelrange(azm, elv, rng, location, pvt.time)

# Handle sequence of PVTs
# tspvcs = [convert._tspvc(p) for p in pvt] # Put into _tspvc

###################################
#### Position from observation ####
###################################

def eci(observation, name="topo loc", forceenv=force.deffe):
    '''Find the ECI position and time of the observations made from
    the location. If observation is a Time or multiple times, find the site
    vector(s). Uses Orekit.'''
    if type(observation) is obs.EarthObservationT:
        tf = _topoframe(observation.loc, name, forceenv)
        obork = units.orig_change_units(observation.obs, units.orkunits)
        tf = TopocentricFrame(forceenv['earth'], \
                              tf.pointAtDistance(
                                  float(obork['azim'].value), \
                                  float(obork['elev'].value), \
                                  float(obork['range'].value)), \
                              name)
        dttm = observation.time
        # The velocity part of this I guess is the frame velocity at
        # that point; of course, it is not the satellite's velocity.
        #
        return convert._pvt(tf, getpvcargs=[dttm, forceenv['celestialframe']], additional=dttm)
        # pos = convert._pvt(tf, getpvcargs=[dttm, forceenv['celestialframe']]).pv['position']
        # return (pos, dttm)


    # This handles times to return the site vector; maybe this should be a different function?
    # elif observation==None:
    #     dttm = atime.abstime(0)
    # else:
    #     dttm = observation

    # If `observation ` is a set of observations
    # if type(dttm.value) is np.ndarray:
    #     arr = [eci(loc, dt) for dt in dttm]
    #     datdict = {name: [ob.cartesian.xyz for ob in arr]}
    #     times = [ob.obstime for ob in arr]
    #     ts = TimeSeries(time=times, data=datdict)
    #     ts[name].info.format = posvel._pos_format
    #     return ts
    # else:
    #     return (pos, dttm) # posvel.makept(pos, dttm)
