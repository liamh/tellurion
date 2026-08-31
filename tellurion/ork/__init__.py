"""Lazy Orekit/JVM initialization and symbol loading for tellurion.ork."""

from importlib import import_module
from threading import Lock

__all__ = ["init_orekit", "ensure_orekit_initialized", "orekit_available"]

_INIT_LOCK = Lock()
_VM_INITIALIZED = False
_INIT_ERROR = None

_OREKIT_SUBMODULES = (
    "convert",
    "dsst",
    "element",
    "force",
    "geog",
    "iod",
    "jacobian",
    "obs",
    "prop",
    "relative",
    "event.eclipse",
    "event.prop",
    "event.util",
    "event.visibility",
)


def init_orekit():
    """Initialize JVM + Orekit data once, explicitly or on first Orekit use."""
    global _VM_INITIALIZED, _INIT_ERROR
    with _INIT_LOCK:
        if _VM_INITIALIZED:
            return
        if _INIT_ERROR is not None:
            raise RuntimeError(
                "Orekit initialization previously failed. "
                "Call `tellurion.ork.init_orekit()` after fixing your Java/Orekit setup."
            ) from _INIT_ERROR
        try:
            import orekit_jpype as orekit

            orekit.initVM()
            from orekit_jpype.pyhelpers import setup_orekit_data

            setup_orekit_data()
            _VM_INITIALIZED = True
        except Exception as exc:
            _INIT_ERROR = exc
            raise RuntimeError(
                "Failed to initialize Orekit JVM/data. "
                "Install/configure Java + orekit-jpype, then call `tellurion.ork.init_orekit()`."
            ) from exc


def ensure_orekit_initialized():
    """Guard used by Orekit-dependent modules before importing Java classes."""
    init_orekit()


def orekit_available():
    """Return True if Orekit initialization succeeds, otherwise False."""
    try:
        init_orekit()
    except RuntimeError:
        return False
    return True


def _import_ork_submodule(name):
    module = import_module(f"tellurion.ork.{name}")
    globals()[name.rsplit(".", 1)[-1]] = module
    return module


def __getattr__(name):
    if name in _OREKIT_SUBMODULES:
        ensure_orekit_initialized()
        return _import_ork_submodule(name)

    ensure_orekit_initialized()
    for module_name in _OREKIT_SUBMODULES:
        module = _import_ork_submodule(module_name)
        if hasattr(module, name):
            value = getattr(module, name)
            globals()[name] = value
            return value
    raise AttributeError(f"module 'tellurion.ork' has no attribute '{name}'")
