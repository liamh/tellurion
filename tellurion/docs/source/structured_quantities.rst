.. _structured_quantities:

Structured Quantities
=====================

Structured quantities are :class:`~astropy.units.Quantity` objects that have multiple named fields,
each with potentially `different units <https://docs.astropy.org/en/stable/units/structured_units.html>`_. They are useful for representing
collections of related physical quantities, such as position vectors,
phase space coordinates, or observational data with multiple measured quantities.

This page demonstrates how to create and work with structured quantities using
the utilities provided in this package, which offer a more convenient interface
than manually constructing structured arrays.

Quick Start
-----------

The simplest way to create a structured quantity is using :func:`tellurion.make_quantity`
with a dictionary::

    >>> from astropy import units as u
    >>> import tellurion as tell
    >>>
    >>> # Create a structured quantity with multiple fields
    >>> data = {'ra': 45.0, 'dec': 30.0}
    >>> coords = tell.make_quantity(data, 'deg')
    >>> coords
    <Quantity (45., 30.) deg>
    >>> coords['ra']
    <Quantity 45. deg>
    >>> coords['dec']
    <Quantity 30. deg>

Basic Usage
-----------

Single Field Structured Quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can create a structured quantity with a single field::

    >>> data = {'distance': 10.0}
    >>> q = tell.make_quantity(data, 'meter')
    >>> q
    <Quantity (10.,) m>

Multiple Fields with Same Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When all fields share the same unit, simply pass a single unit string::

    >>> data = {'x': 1.0, 'y': 2.0, 'z': 3.0}
    >>> position = tell.make_quantity(data, 'm')
    >>> position
    <Quantity (1., 2., 3.) m>

Multiple Fields with Different Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For fields with different units, pass a dictionary mapping field names to units::

    >>> data = {'position': 100.0, 'velocity': 10.0, 'mass': 5.0}
    >>> units = {'position': 'm', 'velocity': 'm/s', 'mass': 'kg'}
    >>> particle = tell.make_quantity(data, units)
    >>> particle
    <Quantity (100., 10., 5.) (m, m / s, kg)>
    >>> particle['position']
    <Quantity 100. m>
    >>> particle['velocity']
    <Quantity 10. m / s>
    >>> particle['mass']
    <Quantity 5. kg>

Working with Arrays
-------------------

Scalar vs. Array Structured Quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, array values create array structured quantities where each element
is a separate row::

    >>> data = {'x': [1, 2, 3], 'y': [4, 5, 6]}
    >>> points = tell.make_quantity(data, 'm')
    >>> points
    <Quantity [(1., 4.), (2., 5.), (3., 6.)] m>
    >>> len(points)
    3
    >>> points[0]  # First point
    <Quantity (1., 4.) m>
    >>> points['x']  # All x coordinates
    <Quantity [1., 2., 3.] m>

Vector Fields as Single Elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To treat an array as a single field value (e.g., a 3D vector), use ``is_scalar=True``::

    >>> data = {'position': [1.0, 2.0, 3.0]}
    >>> vector = tell.make_quantity(data, 'm', is_scalar=True)
    >>> vector.isscalar
    True
    >>> vector['position']
    <Quantity [1., 2., 3.] m>

This is useful for representing quantities where each "row" contains a vector::

    >>> # Multiple 3D positions
    >>> positions_array = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    >>> data = {'position': positions_array}
    >>> positions = tell.make_quantity(data, 'm')
    >>> positions
    <Quantity [([1., 2., 3.],), ([4., 5., 6.],), ([7., 8., 9.],)] m>
    >>> len(positions)  # 3 positions
    3
    >>> positions[0]['position']  # First position vector
    <Quantity [1., 2., 3.] m>
    >>> positions['position'].shape  # (n_positions, 3)
    (3, 3)

Combining Structured Quantities
--------------------------------

Horizontally Stacking (Adding Fields)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`hstack` to combine structured quantities by adding their fields together.
This is like adding columns to a table::

    >>> import tellurion as tell
    >>>
    >>> # Create separate field groups
    >>> positions = tell.make_quantity({'x': [1, 2], 'y': [3, 4]}, 'm')
    >>> velocities = tell.make_quantity({'vx': [10, 20], 'vy': [30, 40]}, 'm/s')
    >>>
    >>> # Combine into phase space
    >>> phase_space = tell.hstack([positions, velocities])
    >>> phase_space
    <Quantity [(1., 3., 10., 30.), (2., 4., 20., 40.)] (m, m, m / s, m / s)>
    >>> phase_space.dtype.names
    ('x', 'y', 'vx', 'vy')

All quantities must have the same length (number of rows)::

    >>> q1 = tell.make_quantity({'a': [1, 2]}, 'm')
    >>> q2 = tell.make_quantity({'b': [3, 4, 5]}, 'm')  # Different length!
    >>> tell.hstack([q1, q2])
    Traceback (most recent call last):
    ...
    ValueError: All structured quantities must have the same length or all be scalar.

Vertically Stacking (Adding Rows)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`vstack` to combine structured quantities by adding rows together.
This is like adding rows to a table::

    >>> import tellurion as tell
    >>>
    >>> # Create separate observations
    >>> obs1 = tell.make_quantity({'ra': 45.0, 'dec': 30.0}, 'deg')
    >>> obs2 = tell.make_quantity({'ra': 120.0, 'dec': -15.0}, 'deg')
    >>> obs3 = tell.make_quantity({'ra': 200.0, 'dec': 60.0}, 'deg')
    >>>
    >>> # Combine into catalog
    >>> catalog = tell.vstack([obs1, obs2, obs3])
    >>> catalog
    <Quantity [(45., 30.), (120., -15.), (200., 60.)] deg>
    >>> len(catalog)
    3
    >>> catalog['ra']
    <Quantity [45., 120., 200.] deg>

All quantities must have the same structure (same field names and units)::

    >>> q1 = tell.make_quantity({'x': 1, 'y': 2}, 'm')
    >>> q2 = tell.make_quantity({'x': 3, 'z': 4}, 'm')  # Different fields!
    >>> tell.vstack([q1, q2])
    Traceback (most recent call last):
    ...
    ValueError: All structured quantities must have the same unit structure.

Unit System Conversion
----------------------

Using a Unit Lookup Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`change_units` function converts quantities to a target unit system
using physical type lookups. This is particularly useful when you want to
ensure consistent units across your calculations::

    >>> import tellurion as tell
    >>>
    >>> # Create quantity in mixed units
    >>> data = {'height': 100, 'width': 50, 'depth': 2.5}
    >>> units = {'height': 'cm', 'width': 'cm', 'depth': 'm'}
    >>> box = tell.make_quantity(data, units)
    >>>
    >>> # Define target unit system (SI base units)
    >>> SI_UNITS = {'length': 'm', 'time': 's', 'mass': 'kg'}
    >>>
    >>> # Convert to SI
    >>> box_si = tell.change_units(box, SI_UNITS)
    >>> box_si
    <Quantity (1., 0.5, 2.5) m>

Using a Master Unit Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A common pattern is to define a master unit dictionary and use subsets of it::

    >>> # Define all possible units in your system
    >>> UNIT_SYSTEM = {
    ...     'x': 'm', 'y': 'm', 'z': 'm',
    ...     'vx': 'm/s', 'vy': 'm/s', 'vz': 'm/s',
    ...     'mass': 'kg',
    ...     'time': 's',
    ...     'energy': 'J'
    ... }
    >>>
    >>> # Create different quantities using subsets
    >>> position = tell.make_quantity({'x': 1, 'y': 2, 'z': 3}, UNIT_SYSTEM)
    >>> velocity = tell.make_quantity({'vx': 10, 'vy': 20, 'vz': 30}, UNIT_SYSTEM)
    >>> particle = tell.make_quantity({'x': 5, 'mass': 2.5, 'energy': 100}, UNIT_SYSTEM)

The unit dictionary can contain more keys than the value dictionary - extra
keys are simply ignored. This allows you to maintain a single comprehensive
unit system definition.

Conversion Between Structured and Unstructured
-----------------------------------------------

Converting to Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^

Use :func:`quantity_to_dict` to extract fields as a dictionary::

    >>> import tellurion as tell
    >>>
    >>> data = {'ra': 45.0, 'dec': 30.0, 'distance': 100.0}
    >>> units = {'ra': 'deg', 'dec': 'deg', 'distance': 'pc'}
    >>> observation = tell.make_quantity(data, units)
    >>>
    >>> obs_dict = tell.quantity_to_dict(observation)
    >>> obs_dict
    {'ra': <Quantity 45. deg>, 'dec': <Quantity 30. deg>, 'distance': <Quantity 100. pc>}
    >>> obs_dict['ra']
    <Quantity 45. deg>

Converting to Array
^^^^^^^^^^^^^^^^^^^

Use :func:`quantity_to_array` to extract numerical values as a numpy array,
discarding unit information::

    >>> import tellurion as tell
    >>>
    >>> data = {'x': [1, 2, 3], 'y': [4, 5, 6]}
    >>> points = tell.make_quantity(data, 'm')
    >>>
    >>> array = tell.quantity_to_array(points)
    >>> array
    array([[1., 4.],
           [2., 5.],
           [3., 6.]])

This is useful when interfacing with libraries that don't understand astropy units::

    >>> import numpy as np
    >>> # After verifying units are correct, extract for computation
    >>> distances = np.linalg.norm(array, axis=1)
    >>> distances
    array([4.123..., 5.385..., 6.708...])

Real-World Examples
-------------------

Astronomical Catalog
^^^^^^^^^^^^^^^^^^^^

Creating a simple astronomical catalog with positions and magnitudes::

    >>> # Define entries
    >>> sources = {
    ...     'ra': [10.68, 83.63, 201.30],      # degrees
    ...     'dec': [41.27, -5.39, -43.02],     # degrees
    ...     'distance': [0.77, 0.41, 4.37],    # parsecs
    ...     'vmag': [0.03, 0.12, -0.27]        # magnitudes (dimensionless)
    ... }
    >>>
    >>> units = {
    ...     'ra': 'deg',
    ...     'dec': 'deg',
    ...     'distance': 'pc',
    ...     'vmag': ''  # dimensionless
    ... }
    >>>
    >>> catalog = tell.make_quantity(sources, units)
    >>> catalog['ra']
    <Quantity [10.68, 83.63, 201.3] deg>
    >>> catalog['distance']
    <Quantity [0.77, 0.41, 4.37] pc>

Particle Simulation
^^^^^^^^^^^^^^^^^^^

Representing particles with positions, velocities, and masses::

    >>> import numpy as np
    >>>
    >>> n_particles = 100
    >>>
    >>> # Initial conditions
    >>> positions = tell.make_quantity({
    ...     'x': np.random.randn(n_particles),
    ...     'y': np.random.randn(n_particles),
    ...     'z': np.random.randn(n_particles)
    ... }, 'm')
    >>>
    >>> velocities = tell.make_quantity({
    ...     'vx': np.random.randn(n_particles),
    ...     'vy': np.random.randn(n_particles),
    ...     'vz': np.random.randn(n_particles)
    ... }, 'm/s')
    >>>
    >>> masses = tell.make_quantity({
    ...     'mass': np.random.uniform(0.1, 10.0, n_particles)
    ... }, 'kg')
    >>>
    >>> # Combine into single phase space
    >>> particles = tell.hstack([positions, velocities, masses])
    >>>
    >>> # Access properties
    >>> particles['mass']  # All masses
    >>> particles[0]  # First particle (all properties)
    >>> particles['x']  # All x positions

Time Series Data
^^^^^^^^^^^^^^^^

Organizing time-series measurements::

    >>> # Measurement times
    >>> times = np.linspace(0, 10, 50)
    >>>
    >>> # Simulated measurements with noise
    >>> temperature = 20 + 5 * np.sin(times) + np.random.randn(50) * 0.5
    >>> pressure = 101.3 + 2 * np.cos(times) + np.random.randn(50) * 0.1
    >>>
    >>> # Create time series
    >>> measurements = tell.make_quantity({
    ...     'time': times,
    ...     'temperature': temperature,
    ...     'pressure': pressure
    ... }, {
    ...     'time': 's',
    ...     'temperature': 'deg_C',
    ...     'pressure': 'kPa'
    ... })
    >>>
    >>> # Plot (after converting to arrays)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(measurements['time'].value,
    ...          measurements['temperature'].value)
    >>> plt.xlabel(f"Time ({measurements['time'].unit})")
    >>> plt.ylabel(f"Temperature ({measurements['temperature'].unit})")

Comparison with Manual Construction
------------------------------------

The traditional way to create structured quantities in AstroPy requires
manually constructing structured arrays::

    >>> # Traditional approach (verbose)
    >>> import numpy as np
    >>> from astropy import units as u
    >>>
    >>> dt = np.dtype([('x', 'f8'), ('y', 'f8')])
    >>> struct_array = np.array([(1.0, 2.0)], dtype=dt)
    >>> traditional = u.Quantity(struct_array, unit=u.StructuredUnit((u.m, u.m)))

Using ``make_quantity`` is more concise and readable::

    >>> # Using make_quantity (concise)
    >>> import tellurion as tell
    >>>
    >>> modern = tell.make_quantity({'x': 1.0, 'y': 2.0}, 'm')

Both produce equivalent results, but ``make_quantity`` is:

- **More readable**: Uses familiar Python dictionaries
- **Less error-prone**: No manual dtype construction
- **More flexible**: Easy to add/remove fields
- **More pythonic**: Natural syntax for named data

Reference/API
-------------

.. autofunction:: tellurion.make_quantity

.. autofunction:: tellurion.change_units

.. autofunction:: tellurion.hstack

.. autofunction:: tellurion.vstack

.. autofunction:: tellurion.quantity_to_dict

.. autofunction:: tellurion.quantity_to_array

See Also
--------

- :ref:`astropy:unit_equivalencies`
- :ref:`astropy:quantity_arithmetic`
- :ref:`astropy:astropy-units-structured-units` (Original AstroPy documentation)
