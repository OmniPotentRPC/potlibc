#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Check that the committed exhaustive inventory matches vendor metadata."""

import importlib.util
import json
import os
import sys
from pathlib import Path


def check_catalog(
    inventory_path: Path,
    vendor_path: Path,
    expected_path: Path,
    plugin_path: Path,
) -> None:
    spec = importlib.util.spec_from_file_location(
        "potlibc_inventory", inventory_path
    )
    if spec is None or spec.loader is None:
        raise SystemExit("cannot import surface inventory")
    inventory_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inventory_module)

    actual = inventory_module.inventory(vendor_path)
    expected = json.loads(expected_path.read_text())
    if actual != expected:
        raise SystemExit(
            "data/surfaces.json is stale; regenerate it with "
            "tools/surface_inventory.py"
        )
    if len({surface["id"] for surface in actual["surfaces"]}) != len(
        actual["surfaces"]
    ):
        raise SystemExit("surface identifiers are not unique")
    for surface in actual["surfaces"]:
        module = plugin_path / f"libpotlibc-{surface['id']}.so"
        if not module.is_file():
            raise SystemExit(f"surface has no built module: {surface['id']}")
    for support in actual["non_surfaces"]:
        if not support["owner_surface_ids"]:
            raise SystemExit(
                f"non-surface file has no owning surface: {support['path']}"
            )


def test_catalog() -> None:
    root = Path(__file__).resolve().parents[1]
    plugin_dir = Path(
        os.environ.get("POTLIBC_PLUGIN_DIR", str(root / "build"))
    ).resolve()
    check_catalog(
        root / "tools" / "surface_inventory.py",
        root / "upstream" / "potlib2023",
        root / "data" / "surfaces.json",
        plugin_dir,
    )


if __name__ == "__main__":
    if len(sys.argv) == 5:
        check_catalog(*map(Path, sys.argv[1:]))
    else:
        test_catalog()
