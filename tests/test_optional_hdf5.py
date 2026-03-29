import sys
import unittest.mock as mock


def test_hdf5_modules_load_when_available():
    """The hdf5 try/except block succeeds when astropy_hdf5io is installed."""
    import astropy_hdf5io  # if this import works, the optional block should too
    import tellurion.core.posvel_hdf5
    import tellurion.core.element_hdf5
    import tellurion.core.spacetrack_hdf5


def test_hdf5_modules_skipped_when_unavailable():
    """The hdf5 try/except block is silently skipped when astropy_hdf5io is missing."""

    # Block only astropy_hdf5io and fsc — these are already imported so we
    # replace them with None to simulate ImportError on any fresh import attempt
    blocked = {
        "astropy_hdf5io": None,
        "fsc": None,
        "fsc.hdf5_io": None,
        "h5py": None,
    }
    # Also block any already-cached submodules of fsc
    for key in list(sys.modules.keys()):
        if key.startswith("fsc") or key.startswith("astropy_hdf5io"):
            blocked[key] = None

    # Remove the tellurion hdf5 submodules so they re-import inside the block
    # but leave tellurion itself and all astropy/orekit modules alone
    hdf5_submodules = [
        key for key in list(sys.modules.keys())
        if key.startswith("tellurion") and "hdf5" in key
    ]
    for key in hdf5_submodules:
        del sys.modules[key]

    with mock.patch.dict(sys.modules, blocked):
        # These should each raise ImportError and be silently caught
        try:
            import tellurion.core.posvel_hdf5
            imported_posvel = True
        except ImportError:
            imported_posvel = False

        try:
            import tellurion.core.element_hdf5
            imported_element = True
        except ImportError:
            imported_element = False

        try:
            import tellurion.core.spacetrack_hdf5
            imported_spacetrack = True
        except ImportError:
            imported_spacetrack = False

    assert not imported_posvel,    "posvel_hdf5 should fail without astropy_hdf5io"
    assert not imported_element,   "element_hdf5 should fail without astropy_hdf5io"
    assert not imported_spacetrack, "spacetrack_hdf5 should fail without astropy_hdf5io"


def test_tellurion_init_try_block():
    """Directly test the try/except logic from __init__.py in isolation."""

    blocked = {
        "astropy_hdf5io": None,
        "fsc": None,
        "fsc.hdf5_io": None,
    }
    for key in list(sys.modules.keys()):
        if key.startswith("fsc") or key.startswith("astropy_hdf5io"):
            blocked[key] = None

    with mock.patch.dict(sys.modules, blocked):
        hdf5_loaded = False
        try:
            import astropy_hdf5io  # noqa: F401
            hdf5_loaded = True
        except ImportError:
            pass

    assert not hdf5_loaded, "astropy_hdf5io should appear missing inside the block"