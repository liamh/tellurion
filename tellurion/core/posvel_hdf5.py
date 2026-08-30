"""
HDF5 serialization support for PositionT and PositionVelocityT classes.

This module extends astropy-hdf5io to support the custom position/velocity classes.
Simply import this module after importing astropy_hdf5io to enable serialization.
"""

from fsc.hdf5_io import subscribe_hdf5

from tellurion.core.posvel import PositionT, PositionVelocityT

##################################################
#### HDF5 Support for PositionT ####
##################################################

def _position_t_to_hdf5(self, hdf5_handle):
    """
    Serialize PositionT to HDF5.

    Parameters
    ----------
    hdf5_handle : HDF5 group
        HDF5 group to write data to
    """
    hdf5_handle["type_tag"] = "tellurion.core.posvel.PositionT"

    # Store the cartesian coordinates
    from fsc.hdf5_io import to_hdf5
    to_hdf5(self.cartesian, hdf5_handle.create_group("cartesian"))

    # Store time if present
    if self.time is not None:
        to_hdf5(self.time, hdf5_handle.create_group("time"))
        hdf5_handle.attrs["has_time"] = True
    else:
        hdf5_handle.attrs["has_time"] = False

    # Store aux dictionary
    if self.aux:
        hdf5_handle.attrs["has_aux"] = True
        # Store aux as attributes (simple key-value pairs)
        aux_group = hdf5_handle.create_group("aux")
        for key, value in self.aux.items():
            # Try to store as HDF5 object, fallback to string
            try:
                to_hdf5(value, aux_group.create_group(str(key)))
            except (TypeError, ValueError):
                aux_group.attrs[str(key)] = str(value)
    else:
        hdf5_handle.attrs["has_aux"] = False

    # Store the isscalar flag
    hdf5_handle.attrs["isscalar"] = self.isscalar


# Monkey-patch the to_hdf5 method onto PositionT
PositionT.to_hdf5 = _position_t_to_hdf5


@subscribe_hdf5("tellurion.core.posvel.PositionT", check_on_load=False)
class _PositionTDeserializer:
    """Deserializer for PositionT"""

    @classmethod  # This is critical!
    def from_hdf5(cls, hdf5_handle):
        """
        Deserialize PositionT from HDF5.

        Parameters
        ----------
        hdf5_handle : HDF5 group
            HDF5 group to read data from

        Returns
        -------
        PositionT
            Reconstructed PositionT object
        """
        from fsc.hdf5_io import from_hdf5

        # Load cartesian coordinates
        cartesian = from_hdf5(hdf5_handle["cartesian"])

        # Load time if present
        if hdf5_handle.attrs["has_time"]:
            time = from_hdf5(hdf5_handle["time"])
        else:
            time = None

        # Load aux dictionary
        aux = {}
        if hdf5_handle.attrs.get("has_aux", False):
            aux_group = hdf5_handle["aux"]
            # Try to load from subgroups first
            for key in aux_group.keys():
                try:
                    aux[key] = from_hdf5(aux_group[key])
                except (TypeError, ValueError, KeyError):
                    pass
            # Also load from attributes
            for key in aux_group.attrs.keys():
                if key not in aux:
                    aux[key] = aux_group.attrs[key]

        # Create the PositionT object
        return PositionT(
            time=time,
            cartesian=cartesian,
            aux=aux
        )


##################################################
#### HDF5 Support for PositionVelocityT ####
##################################################

def _position_velocity_t_to_hdf5(self, hdf5_handle):
    """
    Serialize PositionVelocityT to HDF5.

    Parameters
    ----------
    hdf5_handle : HDF5 group
        HDF5 group to write data to
    """
    hdf5_handle["type_tag"] = "tellurion.core.posvel.PositionVelocityT"

    # Store the cartesian state vector (position + velocity)
    from fsc.hdf5_io import to_hdf5
    to_hdf5(self.cartesian, hdf5_handle.create_group("cartesian"))

    # Store time if present
    if self.time is not None:
        to_hdf5(self.time, hdf5_handle.create_group("time"))
        hdf5_handle.attrs["has_time"] = True
    else:
        hdf5_handle.attrs["has_time"] = False

    # Store aux dictionary
    if self.aux:
        hdf5_handle.attrs["has_aux"] = True
        aux_group = hdf5_handle.create_group("aux")
        for key, value in self.aux.items():
            try:
                to_hdf5(value, aux_group.create_group(str(key)))
            except (TypeError, ValueError):
                aux_group.attrs[str(key)] = str(value)
    else:
        hdf5_handle.attrs["has_aux"] = False

    # Store the isscalar flag
    hdf5_handle.attrs["isscalar"] = self.isscalar


# Monkey-patch the to_hdf5 method onto PositionVelocityT
PositionVelocityT.to_hdf5 = _position_velocity_t_to_hdf5


@subscribe_hdf5("tellurion.core.posvel.PositionVelocityT", check_on_load=False)
class _PositionVelocityTDeserializer:
    """Deserializer for PositionVelocityT"""

    @classmethod  # This is critical!
    def from_hdf5(cls, hdf5_handle):
        """
        Deserialize PositionVelocityT from HDF5.

        Parameters
        ----------
        hdf5_handle : HDF5 group
            HDF5 group to read data from

        Returns
        -------
        PositionVelocityT
            Reconstructed PositionVelocityT object
        """
        from fsc.hdf5_io import from_hdf5

        # Load cartesian state vector
        cartesian = from_hdf5(hdf5_handle["cartesian"])

        # Load time if present
        if hdf5_handle.attrs["has_time"]:
            time = from_hdf5(hdf5_handle["time"])
        else:
            time = None

        # Load aux dictionary
        aux = {}
        if hdf5_handle.attrs.get("has_aux", False):
            aux_group = hdf5_handle["aux"]
            # Try to load from subgroups first
            for key in aux_group.keys():
                try:
                    aux[key] = from_hdf5(aux_group[key])
                except (TypeError, ValueError, KeyError):
                    pass
            # Also load from attributes
            for key in aux_group.attrs.keys():
                if key not in aux:
                    aux[key] = aux_group.attrs[key]

        # Create the PositionVelocityT object
        return PositionVelocityT(
            time=time,
            cartesian=cartesian,
            aux=aux
        )
