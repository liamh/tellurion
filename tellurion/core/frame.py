"""Typed coordinate-frame metadata for Tellurion state objects."""

import dataclasses
from enum import Enum


class FrameKind(str, Enum):
    """Classify a coordinate frame by how its origin and axes are defined."""

    INERTIAL = "inertial"
    BODY_FIXED = "body-fixed"
    LOCAL_ORBITAL = "local-orbital"
    CUSTOM = "custom"


@dataclasses.dataclass(frozen=True, slots=True)
class FrameDescriptor:
    """Describe the origin and axis convention of a state vector.

    This value object is independent of Orekit and is safe to inspect or
    serialize. Orekit ``Frame`` instances should be resolved from this
    metadata at the integration boundary rather than stored here.

    Parameters
    ----------
    origin : str
        Origin of the coordinates, such as ``"Earth"`` or
        ``"Solar System Barycenter"``.
    axes : str
        Axis convention, such as ``"EME2000"``, ``"ICRF"``, ``"LVLH"``,
        or ``"NTW"``.
    kind : FrameKind, optional
        Whether the frame is inertial, body-fixed, local orbital, or custom.
    reference : str, optional
        Identifier for the object defining a local frame, when available.
        This is descriptive metadata; it is not the reference trajectory.
    """

    origin: str
    axes: str
    kind: FrameKind = FrameKind.INERTIAL
    reference: str | None = None

    def __post_init__(self):
        if not isinstance(self.origin, str) or not self.origin.strip():
            raise ValueError("Frame origin must be a non-empty string")
        if not isinstance(self.axes, str) or not self.axes.strip():
            raise ValueError("Frame axes must be a non-empty string")
        if not isinstance(self.kind, FrameKind):
            try:
                object.__setattr__(self, "kind", FrameKind(self.kind))
            except ValueError as exc:
                raise ValueError(f"Unknown frame kind: {self.kind!r}") from exc
        if self.reference is not None and (
            not isinstance(self.reference, str) or not self.reference.strip()
        ):
            raise ValueError("Frame reference must be a non-empty string or None")
