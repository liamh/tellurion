from __future__ import annotations

import os
import shutil
import zipfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve

from platformdirs import user_cache_dir

# You can change this default if you prefer a different source.
# Keep this as a ZIP URL so we can extract it directly.
DEFAULT_OREKIT_DATA_URL = (
    "https://gitlab.orekit.org/orekit/orekit-data/-/archive/main/orekit-data-main.zip"
)

ENV_VAR = "OREKIT_DATA_PATH"
APP_NAME = "tellurion"


def _default_data_dir() -> Path:
    return Path(user_cache_dir(APP_NAME)) / "orekit-data"


def get_orekit_data_path() -> Path | None:
    """
    Return configured Orekit data directory if it exists, else None.
    """
    env = os.getenv(ENV_VAR)
    if env:
        p = Path(env).expanduser().resolve()
        if p.exists() and p.is_dir():
            return p

    p = _default_data_dir()
    if p.exists() and p.is_dir():
        return p

    return None


def ensure_orekit_data(
    *,
    auto_download: bool = True,
    url: str = DEFAULT_OREKIT_DATA_URL,
    force: bool = False,
) -> Path:
    """
    Ensure Orekit data exists locally and return its directory path.

    Behavior:
    - If OREKIT_DATA_PATH points to an existing directory, use it.
    - Else if cached data exists, use it.
    - Else optionally download/unpack data into the user cache.
    """
    existing = get_orekit_data_path()
    if existing and not force:
        return existing

    if not auto_download:
        raise RuntimeError(_missing_data_message())

    cache_root = Path(user_cache_dir(APP_NAME))
    cache_root.mkdir(parents=True, exist_ok=True)

    target = _default_data_dir()
    zip_path = cache_root / "orekit-data.zip"
    tmp_extract = cache_root / "orekit-data-tmp"

    # Clean tmp from prior interrupted runs
    if tmp_extract.exists():
        shutil.rmtree(tmp_extract)

    try:
        print("Downloading Orekit data file, this may take some time...", flush=True)
        urlretrieve(url, zip_path)
    except URLError as e:
        raise RuntimeError(
            f"Could not download Orekit data from {url}.\n{_missing_data_message()}"
        ) from e

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(tmp_extract)

    # GitLab archive usually has a single top-level dir like orekit-data-master/
    extracted_dirs = [p for p in tmp_extract.iterdir() if p.is_dir()]
    if len(extracted_dirs) == 1:
        extracted_root = extracted_dirs[0]
    else:
        extracted_root = tmp_extract

    # Replace target atomically-ish
    if target.exists():
        shutil.rmtree(target)
    shutil.move(str(extracted_root), str(target))

    # Cleanup
    if zip_path.exists():
        zip_path.unlink(missing_ok=True)
    if tmp_extract.exists():
        shutil.rmtree(tmp_extract, ignore_errors=True)

    return target


def orekit_data_setup_instructions() -> str:
    """
    Human-readable setup instructions for users.
    """
    return f"""Orekit data is required.
Options:
1) Recommended: let Tellurion download it automatically:
   >>> import tellurion as tr
   >>> tr.ensure_orekit_data()

2) Manual:
   - Download orekit-data
   - Set environment variable {ENV_VAR} to the extracted directory
"""


def _missing_data_message() -> str:
    return (
        "Orekit data directory not found.\n"
        f"Set {ENV_VAR} to a valid orekit-data directory, or run:\n"
        ">>> import tellurion as tr\n"
        ">>> tr.ensure_orekit_data()\n"
    )
