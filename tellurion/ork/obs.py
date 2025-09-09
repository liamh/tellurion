from org.orekit.frames import TopocentricFrame
from tellurion.core import astro
from tellurion.core import posvel
from tellurion.core import obs
from tellurion.ork import force
from tellurion.ork import convert
from tellurion.ork import geog

######################################
#### Observation values from PVT  ####
######################################

# Make use of element.statefnval
# element.statefnval(pvt, ['azim', 'elev'], obsdict, location)

def _topoframe(loc, name, forceenv=force.deffe):
    return TopocentricFrame(forceenv['earth'], geog.geodpt(loc), name)

# Example demoa.prop.eclipse.ephempvt

def aer(pvt, location, forceenv=force.deffe):
    '''The az-el-range for a PVT'''
    tspvc = convert._tspvc(pvt)
    oktime = tspvc.getDate()
    okpos = tspvc.getPosition()
    tf = _topoframe(location, "local frame", forceenv)
    cf = forceenv['celestialframe']
    azm = tf.getAzimuth(okpos, cf, oktime)*astro.orkunits["angle"]
    elv = tf.getElevation(okpos, cf, oktime)*astro.orkunits["angle"]
    rng = tf.getRange(okpos, cf, oktime)*astro.orkunits["length"]
    return azelrange(azm, elv, rng, location, pvt.time)

# Handle sequence of PVTs
# tspvcs = [convert._tspvc(p) for p in pvt] # Put into _tspvc

##################################
#### Make observation object  ####
##################################
