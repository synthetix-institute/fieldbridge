"""Shared fixtures.

The atlas index is not tracked in git: it is a ~16 MB generated artifact
attached to a release (see docs/NEW_FIELD.md). A fresh clone therefore has no
atlas, and every atlas-dependent test must say so rather than fail, so that a
red suite still means broken code.
"""

from __future__ import annotations

import pytest

from fieldbridge.database import default_data_dir


def atlas_present() -> bool:
    return (default_data_dir() / "index" / "hyperion_static_index.json").exists()


requires_atlas = pytest.mark.skipif(
    not atlas_present(),
    reason="atlas index absent; fetch it with scripts/fetch_atlas.py",
)
