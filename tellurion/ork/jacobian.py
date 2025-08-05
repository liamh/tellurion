import numpy as np
from org.hipparchus.linear import MatrixUtils, RealMatrix
from org.orekit.propagation import MatricesHarvester

# setupMatricesComputation https://www.orekit.org/static/apidocs/org/orekit/propagation/AbstractPropagator.html#setupMatricesComputation(java.lang.String,org.hipparchus.linear.RealMatrix,org.orekit.utils.DoubleArrayDictionary)
# Tutorial example https://gitlab.orekit.org/orekit/orekit-tutorials/-/blob/develop/src/main/java/org/orekit/tutorials/propagation/StateTransitionMatrixPropagation.java

def _add_stm(propagator, dimension):
    '''Request that the state-transition matrix be added to the propagation computation; dimension = 6 without mass, 7 with mass. Returns the harvester.'''
    initstm = MatrixUtils.createRealIdentityMatrix(dimension); # or could be None if dimension = 6
    return propagator.setupMatricesComputation("STM", initstm, None);

def stm(harvester, finalstate):
    '''Return the computed state-transition matrix; note that there are no units conversion, this is in SI units per Orekit standard.'''
    return np.array(harvester.getStateTransitionMatrix(finalstate).getData())
