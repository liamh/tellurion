"""Start the astrodynamics package
        import astrodynamics as ad
        import astrodynamics.example as exmp
        ad.propagate.prop(exmp.ex1,300.0,maximum_tof=86400.0)
        ad.propagate.prop(exmp.ex1,[1200.0,1500.0,1800.0])
   Complete example is exmp.propdemo().
"""
# Define the __all__ variable
# __all__ = ["propagate", "orbit", "frames"]

# Import the submodules
from . import orbit # imports from here: astro, force, ecis, orbit, dttm
from . import propagate # cartprodparam, ecis, orbit
