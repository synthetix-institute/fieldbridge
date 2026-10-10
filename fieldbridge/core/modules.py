"""The module families of FieldBridge and their command-line entry points (shared layer M0).

Each entry names a family, the module that provides ``add_parser(sub)`` for its subcommand, and the installation extra
it needs. ``mount`` adds them to the command line in order; a family whose dependencies are not installed (an
ImportError on import) is left out, as before.
"""
from __future__ import annotations

import importlib
from typing import Tuple

MODULES: Tuple[Tuple[str, str, str], ...] = (
    ("memory", "fieldbridge.memory.cli", "memory"),
    ("quantum", "fieldbridge.quantum.cli", "construction"),
    ("regulation", "fieldbridge.regulation.cli", "memory"),
    ("computation", "fieldbridge.computation.cli", "memory"),
    ("heredity", "fieldbridge.heredity.cli", "memory"),
    ("decision", "fieldbridge.decision.cli", "memory"),
)


def mount(sub) -> list:
    """Add the subcommand of every family that can be imported; returns the names mounted."""
    mounted = []
    for name, path, _extra in MODULES:
        try:
            module = importlib.import_module(path)
        except ImportError:
            continue
        module.add_parser(sub)
        mounted.append(name)
    return mounted
