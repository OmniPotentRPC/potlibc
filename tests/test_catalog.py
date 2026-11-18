#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Check that the committed exhaustive inventory matches vendor metadata."""

import importlib.util
import json
import sys
from pathlib import Path

inventory_path, vendor_path, expected_path, plugin_path = map(Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location("potlibc_inventory", inventory_path)
if spec is None or spec.loader is None:
    raise SystemExit("cannot import surface inventory")
inventory_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory_module)

actual = inventory_module.inventory(vendor_path)
expected = json.loads(expected_path.read_text())


def require_ch52008_gradient_scale(catalog: dict, label: str) -> None:
    matches = [
        surface for surface in catalog["surfaces"] if surface["id"] == "ch52008"
    ]
    if len(matches) != 1:
        raise SystemExit(f"{label} must contain exactly one ch52008 record")
    scale = matches[0].get("gradient_scale")
    if scale != inventory_module.CH52008_GRADIENT_SCALE:
        raise SystemExit(
            f"{label} must record gradient_scale for ch52008 (got {scale!r})"
        )


require_ch52008_gradient_scale(actual, "regenerated inventory")
require_ch52008_gradient_scale(expected, "data/surfaces.json")
if actual != expected:
    raise SystemExit(
        "data/surfaces.json is stale; regenerate it with "
        "tools/surface_inventory.py"
    )
for surface in actual["surfaces"]:
    scale = surface.get("gradient_scale")
    if scale is None:
        continue
    if not isinstance(scale, (int, float)) or scale <= 0:
        raise SystemExit("gradient_scale values must be positive numbers")
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
        raise SystemExit(f"non-surface file has no owning surface: {support['path']}")
