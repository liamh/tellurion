#!/usr/bin/env python3
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path


TUTORIAL_DIR = Path("docs/source/tutorials")
DOCS_BASE = "https://tellurion.readthedocs.io/en/stable"

# Match MyST/Sphinx Python roles in markdown:
# {py:class}`~tellurion.PositionVelocityT`
# {py:func}`tellurion.abstime`
# Support both {py:class}`...` and {{py:class}}`...`
ROLE_RE = re.compile(r"\{\{?py:(class|func|meth|attr|mod)\}?\}`(~?)([^`]+)`")

def role_to_markdown(target: str, shorten: bool) -> str:
    label = target.split(".")[-1] if shorten else target
    href = f"{DOCS_BASE}/api/tellurion.html#{target}"
    return f"[{label}]({href})"


def rewrite_markdown(text: str) -> str:
    def repl(m: re.Match) -> str:
        _kind = m.group(1)
        shorten = m.group(2) == "~"
        target = m.group(3).strip()
        return role_to_markdown(target, shorten)

    return ROLE_RE.sub(repl, text)


def postprocess_ipynb(path: Path) -> None:
    nb = json.loads(path.read_text(encoding="utf-8"))
    changed = False

    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "markdown":
            continue
        src = cell.get("source", [])
        text = "".join(src) if isinstance(src, list) else str(src)
        new_text = rewrite_markdown(text)
        if new_text != text:
            cell["source"] = [new_text]
            changed = True

    if changed:
        path.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        print(f"updated: {path}")
    else:
        print(f"no changes: {path}")


def main() -> None:
    # 1) Convert all percent-format tutorial .py files to .ipynb
    py_files = sorted(TUTORIAL_DIR.glob("*.py"))
    if not py_files:
        raise SystemExit(f"No tutorial .py files found in {TUTORIAL_DIR}")

    subprocess.run(
        ["jupytext", "--to", "ipynb", *[str(p) for p in py_files]],
        check=True,
    )

    # 2) Rewrite roles to plain markdown links in notebooks for Colab
    for ipynb in sorted(TUTORIAL_DIR.glob("*.ipynb")):
        postprocess_ipynb(ipynb)


if __name__ == "__main__":
    main()
