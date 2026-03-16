"""
Script to create pre-recorded HDF5 fixtures for spacetrack_hdf5 tests.

Run this script ONCE when you have access to space-track.org:

    python tests/make_spacetrack_fixtures.py

It saves a dict of MeanElementSetT objects to tests/data/spacetrack_fixtures.h5
so that TestSpacetrackHdf5Fixture in test_spacetrack_hdf5.py can run offline.

Requires stclient to be configured in your IPython startup script:
    # ~/.ipython/profile_default/startup/50-spacetrack.py
    import spacetrack
    stclient = spacetrack.SpaceTrackClient(identity="myemail@example.com",
                                           password="mypw")

Or set it up directly here before running.
"""
import pathlib
import astropy_hdf5io       # registers astropy serializers
from fsc.hdf5_io import save
import tellurion as tell

DATA_DIR = pathlib.Path(__file__).parent / 'data'
DATA_DIR.mkdir(exist_ok=True)
FIXTURE_FILE = DATA_DIR / 'spacetrack_fixtures.h5'

SATNUMS = [41335, 25544]   # Sentinel-3A, ISS

# Uncomment and fill in if not using the IPython startup script:
# stclient = st.SpaceTrackClient(identity="myemail@example.com", password="mypw")

elements = tell.spacetrack_latest(stclient, SATNUMS)  # noqa: F821  (defined in startup)
save(elements, str(FIXTURE_FILE))

print(f"Saved fixtures to {FIXTURE_FILE}")
for name, obj in elements.items():
    print(f"  {name}: catid={obj.scdata['catid']}, epoch={obj.t.iso}")
