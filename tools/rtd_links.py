#!/usr/bin/env python3
"""Replace MyST/Sphinx Python roles in notebooks with Read the Docs links.

Usage: python tools/rtd_links.py docs/source/tutorials/*.ipynb
"""
import re
import sys
import time
import urllib.error
import urllib.request
import zlib

DOCS_BASE = "https://tellurion.readthedocs.io/en/latest"
ROLE_RE = re.compile(r"\{(?:py:)?(class|meth|func|attr|mod|exc|data)\}`(~?)([^`]+)`")
INV_LINE = re.compile(r"(.+?)\s+(\S+)\s+(-?\d+)\s+(\S*)\s+(.*)")
HEADERS = {
    # The default Python-urllib User-Agent is often rejected with 403.
    "User-Agent": "Mozilla/5.0 (compatible; tellurion-docs-linker/1.0; "
                  "+https://github.com/liamh/tellurion)",
    "Accept": "*/*",
}


def fetch(url, attempts=4):
    """GET url with a browser-like User-Agent, retrying on transient errors."""
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            return urllib.request.urlopen(req, timeout=30).read()
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            print(f"fetch attempt {i + 1}/{attempts} failed: {exc}",
                  file=sys.stderr)
            time.sleep(2 * (i + 1))
    raise SystemExit(f"Could not download {url}: {last}")


def load_inventory(base):
    """Return {object name: absolute URL} for Python-domain objects."""
    raw = fetch(f"{base}/objects.inv")
    # Header is 4 text lines, then zlib-compressed entries.
    pos = 0
    for _ in range(4):
        pos = raw.index(b"\n", pos) + 1
    entries = {}
    for line in zlib.decompress(raw[pos:]).decode("utf-8").splitlines():
        m = INV_LINE.match(line)
        if not m:
            continue
        name, domain_role, _prio, uri, _disp = m.groups()
        if not domain_role.startswith("py:"):
            continue
        entries.setdefault(name, f"{base}/{uri.replace('$', name)}")
    return entries


def make_replacer(inv, missing):
    def repl(m):
        kind, tilde, target = m.group(1), m.group(2), m.group(3).strip()
        label = target.split(".")[-1] if tilde else target
        if kind in ("meth", "func"):
            label += "()"
        url = inv.get(target)
        if url is None:
            missing.add(target)
            return f"`{label}`"  # plain code text, never a broken link
        return f"[{label}]({url})"
    return repl


def main(paths):
    inv = load_inventory(DOCS_BASE)
    print(f"Loaded {len(inv)} Python objects from {DOCS_BASE}/objects.inv")
    missing = set()
    for path in paths:
        text = open(path, encoding="utf-8").read()
        new = ROLE_RE.sub(make_replacer(inv, missing), text)
        if new != text:
            open(path, "w", encoding="utf-8").write(new)
            print(f"updated: {path}")
    if missing:
        print("Not found in objects.inv (left as plain text):", file=sys.stderr)
        for t in sorted(missing):
            print(f"  {t}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
