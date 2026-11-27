import sys

import numpy as np
import pytest
import astropy.units as u

import tellurion as tell
import tellurion.core.oem as oem_reader


def _segment(center, axes, start, stop, rows):
    states = "\n".join(
        f"{epoch} {' '.join(str(value) for value in values)}"
        for epoch, values in rows
    )
    return f"""META_START
OBJECT_NAME = TEST SATELLITE
OBJECT_ID = 2025-001A
CENTER_NAME = {center}
REF_FRAME = {axes}
TIME_SYSTEM = UTC
START_TIME = {start}
STOP_TIME = {stop}
META_STOP
{states}
"""


def _write_oem(path, segments):
    path.write_text(
        "CCSDS_OEM_VERS = 2.0\n"
        "CREATION_DATE = 2025-01-01T00:00:00\n"
        "ORIGINATOR = Tellurion tests\n\n"
        + "\n".join(segments),
        encoding="ascii",
    )
    return path


def _two_state_rows(start, stop):
    return [
        (start, [1, 2, 3, 4, 5, 6]),
        (stop, [2, 3, 4, 5, 6, 7]),
    ]


def test_oem_single_segment_has_frame_and_standard_units(tmp_path):
    pytest.importorskip("oem")
    path = _write_oem(
        tmp_path / "earth.oem",
        [
            _segment(
                "EARTH", "EME2000", "2025-01-01T00:00:00",
                "2025-01-01T00:00:10",
                _two_state_rows("2025-01-01T00:00:00", "2025-01-01T00:00:10"),
            )
        ],
    )

    pvt = tell.oempvt(path)

    assert isinstance(pvt, tell.PositionVelocityT)
    assert pvt.frame == tell.FrameDescriptor("Earth", "EME2000")
    assert len(pvt) == 2
    np.testing.assert_allclose(pvt.position_vector[0].to_value(u.km), [1, 2, 3])
    np.testing.assert_allclose(
        pvt.velocity_vector[0].to_value(u.km / u.s), [4, 5, 6]
    )


def test_oem_frame_changes_return_tuple_and_preserve_unknowns(tmp_path):
    pytest.importorskip("oem")
    first = ("2025-01-01T00:00:00", "2025-01-01T00:00:10")
    second = ("2025-01-01T00:00:20", "2025-01-01T00:00:30")
    third = ("2025-01-01T00:00:40", "2025-01-01T00:00:50")
    path = _write_oem(
        tmp_path / "mixed.oem",
        [
            _segment("EARTH", "EME2000", *first, _two_state_rows(*first)),
            _segment("MARS", "UNKNOWN_AXES", *second, _two_state_rows(*second)),
            _segment("MARS", "UNKNOWN_AXES", *third, _two_state_rows(*third)),
        ],
    )

    result = tell.oempvt(path)

    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result[0].frame == tell.FrameDescriptor("Earth", "EME2000")
    assert result[1].frame == tell.FrameDescriptor(
        "MARS", "UNKNOWN_AXES", kind=tell.FrameKind.CUSTOM
    )
    assert len(result[0]) == 2
    assert len(result[1]) == 4
    assert result[1].time[0] < result[1].time[-1]


def test_oem_infers_body_fixed_and_local_orbital_frames(tmp_path):
    pytest.importorskip("oem")
    start, stop = "2025-01-01T00:00:00", "2025-01-01T00:00:10"
    for axes, kind in (
        ("ITRF-2014", tell.FrameKind.BODY_FIXED),
        ("TNW", tell.FrameKind.LOCAL_ORBITAL),
    ):
        path = _write_oem(
            tmp_path / f"{axes}.oem",
            [
                _segment(
                    "LUNA", axes, start, stop, _two_state_rows(start, stop)
                )
            ],
        )
        pvt = tell.oempvt(path)
        assert pvt.frame == tell.FrameDescriptor("Moon", axes, kind=kind)


def test_oem_dependency_is_loaded_only_when_reader_is_called(monkeypatch, tmp_path):
    pytest.importorskip("oem")
    original = oem_reader.importlib.import_module

    def missing_oem(name):
        if name == "oem":
            raise ModuleNotFoundError("missing", name="oem")
        return original(name)

    monkeypatch.setattr(oem_reader.importlib, "import_module", missing_oem)
    with pytest.raises(ImportError, match=r"tellurion\[oem\]"):
        tell.oempvt(tmp_path / "unused.oem")


def test_tellurion_import_does_not_import_optional_oem():
    code = "import sys, tellurion; assert 'oem' not in sys.modules"
    import subprocess

    subprocess.run([sys.executable, "-c", code], check=True)
