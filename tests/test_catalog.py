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
gradient_scales = {
    surface["id"]: surface.pop("gradient_scale")
    for surface in expected["surfaces"]
    if "gradient_scale" in surface
}
if actual != expected:
    raise SystemExit(
        "data/surfaces.json is stale; regenerate it with "
        "tools/surface_inventory.py"
    )
if any(
    not isinstance(scale, (int, float)) or scale <= 0
    for scale in gradient_scales.values()
):
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
