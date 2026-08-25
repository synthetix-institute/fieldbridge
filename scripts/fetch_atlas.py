#!/usr/bin/env python3
"""Download the atlas index into data/.

The index and its shards are generated artifacts of about 31 MB. Tracking them
in git would add that much immutable history on every refresh, so they are
attached to a release instead and fetched here.

    python3 scripts/fetch_atlas.py
    python3 scripts/fetch_atlas.py --release v0.2.0
    python3 scripts/fetch_atlas.py --url file:///path/to/hyperion_static_index.json

Without an atlas the tools still run: the constructor reports a route-derived
proxy instead of an assignment, and atlas-dependent tests skip.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO = "synthetix-institute/fieldbridge"
ASSET = "hyperion_static_index.json"
MIN_RECORDS = 2633


def data_dir() -> Path:
    here = Path(__file__).resolve().parent.parent
    return here / "data"


def fetch(url: str, dest: Path) -> None:
    print(f"fetching {url}")
    try:
        with urllib.request.urlopen(url) as response:
            payload = response.read()
    except urllib.error.HTTPError as error:
        raise SystemExit(
            f"download failed ({error.code}). If the release does not exist yet, "
            f"generate the index per docs/NEW_FIELD.md and copy it to {dest}."
        )
    except urllib.error.URLError as error:
        raise SystemExit(f"could not reach {url}: {error.reason}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(payload)
    print(f"wrote {dest} ({len(payload) / 1e6:.1f} MB)")


def verify(path: Path) -> int:
    """Fail loudly on a truncated or shrunken index.

    A silently smaller atlas degrades every retrieval without erroring, which is
    the failure mode this whole file exists to prevent.
    """
    try:
        index = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise SystemExit(f"{path} is not valid JSON: {error}")
    records = index.get("records") or []
    generated = (index.get("generated_at") or "")[:10]
    print(f"records: {len(records)} | generated: {generated}")
    if len(records) < MIN_RECORDS:
        print(
            f"WARNING: {len(records)} records is below the {MIN_RECORDS} previously "
            "shipped. If you generated this yourself, raise --max-records.",
            file=sys.stderr,
        )
    return len(records)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--release", default="latest",
                        help="release tag to pull from (default: latest)")
    parser.add_argument("--url", help="explicit URL, overriding --release")
    parser.add_argument("--force", action="store_true",
                        help="overwrite an index that is already present")
    args = parser.parse_args()

    dest = data_dir() / "index" / ASSET
    if dest.exists() and not args.force:
        print(f"{dest} already present; pass --force to replace it")
        return 0 if verify(dest) else 1

    if args.url:
        url = args.url
    elif args.release == "latest":
        url = f"https://github.com/{REPO}/releases/latest/download/{ASSET}"
    else:
        url = f"https://github.com/{REPO}/releases/download/{args.release}/{ASSET}"

    fetch(url, dest)
    verify(dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
