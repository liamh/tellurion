# If astropy-hdf5io is present, all tests will run; if not, the tests
# involving hdf5 will be skipped.

# Run ALL tests (including HDF5)
pytest tests/

# Run with HDF5 mock-absent (skips all hdf5-marked tests)
pytest tests/ -m "not hdf5"

# Run ONLY HDF5 tests
pytest tests/ -m hdf5
