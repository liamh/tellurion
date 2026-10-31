"""Serialization support for tellurion.core.spacetrack.MeanElementSetT

MeanElementSetT is a plain class (not a namedtuple), so fsc.hdf5-io's tuple
interceptor does not apply. The monkey-patched to_hdf5 method is found by
fsc.hdf5-io's to_hdf5() dispatcher, and save_recursive/_is_dict_like correctly
treats it as a serializable leaf (has to_hdf5, not traversed as dict-like).
"""

from fsc.hdf5_io import from_hdf5, subscribe_hdf5, to_hdf5

from tellurion.core.spacetrack import MeanElementSetT

_TYPE_TAG = "tellurion.MeanElementSetT"


def _mean_element_set_to_hdf5(self, hdf5_handle):
    """Serialize MeanElementSetT to HDF5."""
    hdf5_handle["type_tag"] = _TYPE_TAG

    # els: dict of orbital elements (astropy Quantity values)
    els_group = hdf5_handle.create_group("els")
    for key, value in self.els.items():
        to_hdf5(value, els_group.create_group(key))

    # time: astropy.time.Time epoch (scalar or vector)
    to_hdf5(self.time, hdf5_handle.create_group("time"))

    # tle: two ASCII TLE line strings, or two lists of them for a series
    for name, lines in zip(("tle_line1", "tle_line2"), self.tle):
        hdf5_handle.attrs[name] = (lines if isinstance(lines, str)
                                   else [str(x) for x in lines])

    # model: plain string (e.g. 'SGP4')
    hdf5_handle.attrs["model"] = self.model

    # scdata: dict of plain scalars (str / int)
    scdata_group = hdf5_handle.create_group("scdata")
    for key, value in self.scdata.items():
        scdata_group.attrs[key] = value


MeanElementSetT.to_hdf5 = _mean_element_set_to_hdf5


@subscribe_hdf5(_TYPE_TAG, check_on_load=False)
class _MeanElementSetDeserializer:
    """Deserializer for MeanElementSetT."""

    @classmethod
    def from_hdf5(cls, hdf5_handle):
        els = {key: from_hdf5(hdf5_handle["els"][key])
               for key in hdf5_handle["els"].keys()}

        # Files written before the rename stored the epoch as "t"
        time = from_hdf5(hdf5_handle["time" if "time" in hdf5_handle else "t"])

        def _lines(name):
            val = hdf5_handle.attrs[name]
            if isinstance(val, (str, bytes)):
                return str(val)
            return [x.decode() if isinstance(x, bytes) else str(x) for x in val]

        tle = (_lines("tle_line1"), _lines("tle_line2"))

        model = str(hdf5_handle.attrs["model"])

        scdata = {key: hdf5_handle["scdata"].attrs[key]
                  for key in hdf5_handle["scdata"].attrs.keys()}
        if "catid" in scdata:
            scdata["catid"] = int(scdata["catid"])

        return MeanElementSetT(els=els, time=time, tle=tle, model=model, scdata=scdata)