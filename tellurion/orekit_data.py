"""Locate Orekit data: bundled snapshot first, cached download as fallback."""

import os
import shutil
import tempfile
import urllib.error
import urllib.request
import zipfile
from importlib import resources
from pathlib import Path

DATA_URLS = tuple(
    f"https://gitlab.orekit.org/orekit/orekit-data/-/archive/{b}/"
    f"orekit-data-{b}.zip"
    for b in ("main", "master")
)
DATA_URL = DATA_URLS[0]
ENV_VAR = "TELLURION_OREKIT_DATA"
_MARKERS = ("tai-utc.dat", "UTC-TAI.history")


def _valid(path):
    return path.is_dir() and any(any(path.rglob(m)) for m in _MARKERS)


def bundled_data_path():
    """Return the bundled data directory, or None if it is not present."""
    try:
        path = Path(str(resources.files("tellurion") / "_orekit_data"))
    except (ModuleNotFoundError, TypeError):
        return None
    return path if _valid(path) else None


def cache_dir():
    """Directory where downloaded data is cached."""
    env = os.environ.get("TELLURION_CACHE")
    if env:
        return Path(env)
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "tellurion"


def download_orekit_data(dest=None, url=None):
    """Download and extract Orekit data into ``dest``; return its path."""
    dest = Path(dest) if dest else cache_dir() / "orekit-data"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        zpath = Path(tmp) / "orekit-data.zip"
        last_exc = None
        for candidate in (url,) if url else DATA_URLS:
            try:
                with urllib.request.urlopen(candidate, timeout=60) as resp:
                    with open(zpath, "wb") as f:
                        shutil.copyfileobj(resp, f)
                break
            except urllib.error.URLError as exc:
                last_exc = exc
        else:
            raise RuntimeError(
                f"Could not download Orekit data from {DATA_URLS}"
            ) from last_exc
        with zipfile.ZipFile(zpath) as z:
            z.extractall(tmp)
        roots = [p for p in Path(tmp).iterdir() if p.is_dir()]
        if len(roots) != 1 or not _valid(roots[0]):
            raise RuntimeError("Downloaded archive does not contain Orekit data")
        if dest.exists():
            shutil.rmtree(dest)
        shutil.move(str(roots[0]), str(dest))
    return dest


def orekit_data_path():
    """Return a directory with Orekit data.

    Order: ``TELLURION_OREKIT_DATA`` env var, bundled snapshot, cached
    download (fetched if missing).
    """
    env = os.environ.get(ENV_VAR)
    if env:
        return Path(env)
    bundled = bundled_data_path()
    if bundled is not None:
        return bundled
    cached = cache_dir() / "orekit-data"
    if _valid(cached):
        return cached
    return download_orekit_data(cached)
