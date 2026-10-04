"""Lazy Orekit/JVM initialization and symbol loading for tellurion.ork."""

import ast
from importlib import import_module, util as importlib_util
from pathlib import Path
from threading import Lock

__all__ = ["init_orekit", "ensure_orekit_initialized", "orekit_available"]

_INIT_LOCK = Lock()
_SYMBOL_LOCK = Lock()
_VM_INITIALIZED = False
_INIT_ERROR = None

_OREKIT_TOPLEVEL_SUBMODULES = (
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
)
_OREKIT_SCAN_MODULES = _OREKIT_TOPLEVEL_SUBMODULES + (
    "event.eclipse",
    "event.prop",
    "event.util",
    "event.visibility",
)
_OREKIT_SYMBOL_TO_MODULES = None


def init_orekit():
    """Initialize JVM + Orekit data once, explicitly or on first Orekit use."""
    global _VM_INITIALIZED, _INIT_ERROR
    if _VM_INITIALIZED:
        return
    with _INIT_LOCK:
        if _VM_INITIALIZED:
            return
        if _INIT_ERROR is not None:
            raise RuntimeError(
                "Orekit initialization previously failed. "
                "Restart your Python session after fixing your Java/Orekit setup."
            ) from _INIT_ERROR
        try:
            import orekit_jpype as orekit

            orekit.initVM()
            from java.io import File
            from org.orekit.data import DataContext, DirectoryCrawler

            from tellurion.orekit_data import orekit_data_path

            manager = DataContext.getDefault().getDataProvidersManager()
            manager.addProvider(DirectoryCrawler(File(str(orekit_data_path()))))
            _VM_INITIALIZED = True
        except Exception as exc:
            _INIT_ERROR = exc
            raise RuntimeError(
                "Failed to initialize Orekit JVM/data. "
                "Install/configure Java + orekit-jpype, then call "
                "`tellurion.ork.init_orekit()`."
            ) from exc


def ensure_orekit_initialized():
    """Guard used by Orekit-dependent modules before importing Java classes."""
    init_orekit()


def orekit_available():
    """Return True if orekit-jpype bridge module is importable."""
    return importlib_util.find_spec("orekit_jpype") is not None


def _import_ork_submodule(name):
    module = import_module(f"tellurion.ork.{name}")
    if "." not in name:
        globals()[name] = module
    return module


def _discover_exportable_symbols():
    exportable = {}
    base_dir = Path(__file__).resolve().parent
    for module_name in _OREKIT_SCAN_MODULES:
        module_path = base_dir / (module_name.replace(".", "/") + ".py")
        if not module_path.exists():
            continue
        tree = ast.parse(module_path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                if not node.name.startswith("_"):
                    exportable.setdefault(node.name, []).append(module_name)
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and not target.id.startswith("_"):
                        exportable.setdefault(target.id, []).append(module_name)
            elif isinstance(node, ast.AnnAssign):
                target = node.target
                if isinstance(target, ast.Name) and not target.id.startswith("_"):
                    exportable.setdefault(target.id, []).append(module_name)
    return exportable


def __getattr__(name):
    global _OREKIT_SYMBOL_TO_MODULES

    if name in _OREKIT_TOPLEVEL_SUBMODULES:
        return _import_ork_submodule(name)

    if _OREKIT_SYMBOL_TO_MODULES is None:
        with _SYMBOL_LOCK:
            if _OREKIT_SYMBOL_TO_MODULES is None:
                _OREKIT_SYMBOL_TO_MODULES = _discover_exportable_symbols()

    module_names = _OREKIT_SYMBOL_TO_MODULES.get(name)
    if module_names is None:
        raise AttributeError(f"module 'tellurion.ork' has no attribute '{name}'")

    ensure_orekit_initialized()
    for module_name in module_names:
        module = _import_ork_submodule(module_name)
        try:
            value = getattr(module, name)
        except AttributeError:
            continue
        globals()[name] = value
        return value
    raise AttributeError(f"module 'tellurion.ork' has no attribute '{name}'")
