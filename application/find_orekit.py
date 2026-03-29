import types
import jpype
from collections.abc import Mapping, Sequence

def _is_orekit_object(obj):
    """Return True if obj is a JPype/Orekit Java object."""
    try:
        return isinstance(obj, jpype.JClass('java.lang.Object'))
    except Exception:
        return False

def _is_unsaveable(obj):
    """Return True if obj cannot be serialised to HDF5."""
    return (
        _is_orekit_object(obj)
        or isinstance(obj, (
            types.FunctionType,      # regular functions and lambdas
            types.MethodType,        # bound methods
            types.BuiltinFunctionType,
            types.BuiltinMethodType,
        ))
    )

def strip_orekit(obj):
    """
    Recursively walk a munch/dict/list/tuple and remove any objects that
    cannot be serialised to HDF5 (Orekit/JPype objects and functions/lambdas),
    replacing them with None.
    """
    if _is_unsaveable(obj):
        return None

    if isinstance(obj, Mapping):
        cleaned = {k: strip_orekit(v) for k, v in obj.items()}
        try:
            return type(obj)(cleaned)
        except Exception:
            return cleaned

    if type(obj) is list:
        return [strip_orekit(v) for v in obj]

    if type(obj) is tuple:
        return tuple(strip_orekit(v) for v in obj)

    return obj

def find_orekit(obj, path=''):
    """Print the key path of any Orekit objects found in a munch/dict."""
    if _is_orekit_object(obj):
        print(f"  Orekit object at: {path!r}  (type: {type(obj).__name__})")
        return

    if isinstance(obj, Mapping):
        for k, v in obj.items():
            find_orekit(v, f"{path}.{k}" if path else str(k))

    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            find_orekit(v, f"{path}[{i}]")

def find_unsaveable(obj, path=''):
    """Print the key path of any unsaveable objects in a munch/dict."""
    if _is_unsaveable(obj):
        print(f"  Unsaveable at: {path!r}  (type: {type(obj).__name__}, value: {obj!r})")
        return

    if isinstance(obj, Mapping):
        for k, v in obj.items():
            find_unsaveable(v, f"{path}.{k}" if path else str(k))

    elif type(obj) in (list, tuple):
        for i, v in enumerate(obj):
            find_unsaveable(v, f"{path}[{i}]")

# find_unsaveable(demoa)
# find_orekit(demoa)
