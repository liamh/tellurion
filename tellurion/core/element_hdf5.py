"""
HDF5 serialization support for ElementSetT.

This module extends astropy-hdf5io to support the ElementSetT class from
tellurion.core.element. Simply import this module after importing
astropy_hdf5io to enable serialization.

Example
-------
>>> import astropy_hdf5io                   # registers astropy serializers
>>> import tellurion.core.element_hdf5      # registers ElementSetT serializer
>>> from fsc.hdf5_io import save, load
>>> save(my_element_set, 'orbit.hdf5')
>>> loaded = load('orbit.hdf5')
"""

import astropy.units as u
from fsc.hdf5_io import from_hdf5, subscribe_hdf5, to_hdf5

from tellurion.core.element import ElementSetT

##################################################
#### HDF5 Support for ElementSetT            ####
##################################################

def _element_set_t_to_hdf5(self, hdf5_handle):
    """
    Serialize ElementSetT to HDF5.

    Each orbital element field is stored as a raw scalar value plus a unit
    string. This avoids delegating to the astropy-hdf5io Quantity serializer,
    which cannot round-trip 0-d structured-array fields reliably.

    Parameters
    ----------
    hdf5_handle : h5py.Group
        HDF5 group to write data into.
    """
    hdf5_handle["type_tag"] = "tellurion.core.element.ElementSetT"

    # Serialize the epoch Time using the astropy-hdf5io Time serializer
    to_hdf5(self.time, hdf5_handle.create_group("t"))

    # Serialize the structured Quantity field-by-field
    els_grp = hdf5_handle.create_group("els")
    field_names = list(self.elements.dtype.names)
    # Store field order so the structured dtype is rebuilt identically on load
    els_grp.attrs["field_names"] = field_names
    for name in field_names:
        fld_grp = els_grp.create_group(name)
        fld = self.elements[name]
        # Extract the plain Python/numpy scalar and the unit string
        fld_grp["value"] = float(fld.value)
        fld_grp.attrs["unit"] = str(fld.unit)


# Monkey-patch the to_hdf5 method onto ElementSetT
ElementSetT.to_hdf5 = _element_set_t_to_hdf5


@subscribe_hdf5("tellurion.core.element.ElementSetT", check_on_load=False)
class _ElementSetTDeserializer:
    """Deserializer for ElementSetT."""

    @classmethod
    def from_hdf5(cls, hdf5_handle):
        """
        Deserialize ElementSetT from HDF5.

        Parameters
        ----------
        hdf5_handle : h5py.Group
            HDF5 group to read data from.

        Returns
        -------
        ElementSetT
            Reconstructed ElementSetT object.
        """
        # Reconstruct the epoch Time
        t = from_hdf5(hdf5_handle["t"])

        # Reconstruct the structured Quantity field-by-field
        els_grp = hdf5_handle["els"]
        field_names = list(els_grp.attrs["field_names"])

        # Build a dict of plain Quantity objects, then assemble the structured array
        fields = {}
        for name in field_names:
            fld_grp = els_grp[name]
            value = float(fld_grp["value"][()])
            unit  = u.Unit(fld_grp.attrs["unit"])
            fields[name] = u.Quantity(value, unit)

        # Reassemble into a structured Quantity the same way kepler() does,
        # using make_quantity so that units and angle normalisation are applied.
        from tellurion.astro import quantity_utils as quant, units as tunits

        # keppt maps each field name to its physical type.  We infer the
        # phystype from the reconstructed unit so no extra lookup table is needed.
        keppt = {
            "inc": "angle", "argper": "angle", "raan": "angle",
            "ma":  "angle", "ta":     "angle",
            "ecc": "dimensionless", "sma": "length", "memo": "angular speed",
            "radper": "length", "radapo": "length",
            "altper": "length", "altapo": "length",
        }
        phystype_map = {name: keppt[name] for name in field_names if name in keppt}

        isscalar = t.isscalar
        els = quant.make_quantity(fields, phystype_map, isscalar, tunits.prefunits)

        return ElementSetT(els, t)
