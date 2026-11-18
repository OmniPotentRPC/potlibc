#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Ensure triatomic pair-common wrappers follow each surface header order.

References (POTLIB headers):
- FH2 5SEC: G. C. Lynch et al., J. Chem. Phys. 94, 7136-7149 (1991)
- BrH2 LEPS: R. Weston, J. Chem. Phys. 31, 892 (1959)
- ClH2 G3: Allison et al., J. Phys. Chem. (submitted; cited in file header)
"""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "tools"))
from generate_surface import _save_entry_locals  # noqa: E402

surface_expectations = {
    # Header definitions:
    #   R(1) = R(F-H),   R(2) = R(H-H),   R(3) = R(H-F)
    "fh25sec": ("(1,:)", "(2,:)", "(3,:)"),
    #   R(1) = R(Br-H),  R(2) = R(H-H),   R(3) = R(Br-H)
    "brh2leps": ("(1,:)", "(2,:)", "(3,:)"),
    #   R(1) = R(H-Cl),  R(2) = R(H-H'),  R(3) = R(H'-Cl)
    "clh2g3": ("(1,:)", "(2,:)", "(3,:)"),
    "clh2g3v": ("(1,:)", "(2,:)", "(3,:)"),
    #   R(1) = R(Cl-first H), R(2) = R(first H-second H), R(3) = R(Cl-second H)
    "clh2alab": ("(1,:)", "(2,:)", "(3,:)"),
    "hi2kka": ("(1,:)", "(1,:)", "(2,:)"),
    "hf2jotii": ("(1,:)", "(1,:)", "(2,:)"),
    "lifhj": ("(1,:)", "(2,:)", "(1,:)"),
    "lifhjs": ("(1,:)", "(2,:)", "(1,:)"),
    "lifhd": ("(1,:)", "(2,:)", "(1,:)"),
    "nafh1v": ("(1,:)", "(2,:)", "(1,:)"),
}


def test_internal_coordinate_order() -> None:
    build_root = Path(
        os.environ.get("POTLIBC_PLUGIN_DIR", str(root / "build"))
    ).resolve()
    expected = {
        "hi2kka": [("1", "1", "2"), ("2", "1", "3"), ("3", "2", "3")],
        "hf2jotii": [("1", "1", "3"), ("2", "1", "2"), ("3", "2", "3")],
    }
    standard = [("1", "1", "2"), ("2", "2", "3"), ("3", "1", "3")]
    for surface_id in surface_expectations:
        adapter = (build_root / f"{surface_id}-adapter.f90").read_text(
            encoding="utf-8"
        )
        assignments = re.findall(
            r"r\((?:1,)?(\d)\)=distance3\(coords\((\d),:\),coords\((\d),:\)\)",
            adapter,
        )
        if assignments[:3] != expected.get(surface_id, standard):
            raise SystemExit(
                f"{surface_id}: unexpected R order {assignments[:3]!r}"
            )


def test_save_entry_locals_function_headers() -> None:
    headers = [
        "recursive real*4 function f(x)",
        "real function f(x)",
        "character*12 function f(x)",
        "complex function f(x)",
    ]
    for header in headers:
        transformed = _save_entry_locals(
            f"{header}\nentry g(x)\nf=x\nend function f\n"
        )
        if transformed.splitlines()[:2] != [header, "save"]:
            raise SystemExit(f"SAVE inserted outside {header!r}")


def test_all_generated_entry_units_have_save() -> None:
    build_root = Path(
        os.environ.get("POTLIBC_PLUGIN_DIR", str(root / "build"))
    ).resolve()
    catalog = json.loads((root / "data" / "surfaces.json").read_text())
    sources = {
        path.name.removesuffix("-surface.f90"): path
        for path in build_root.glob("*-surface.f90")
    }
    expected_ids = {surface["id"] for surface in catalog["surfaces"]}
    if set(sources) != expected_ids:
        missing = sorted(expected_ids - set(sources))
        raise SystemExit(f"generated sources missing from corpus check: {missing}")
    unit_pattern = re.compile(
        r"(?ims)(^\s*(?:(?:recursive|pure|elemental|impure|module)\s+)*"
        r"(?:(?:double\s+precision|real(?:\s*\*\s*\d+)?|integer|logical|"
        r"character(?:\s*\*\s*\d+)?|complex(?:\s*\*\s*\d+)?|"
        r"type\s*\(\s*dual2\s*\))\s+)?"
        r"(?:subroutine|function)\b.*?"
        r"^\s*end(?:\s+(?:subroutine|function)(?:\s+\w+)?)?\s*$)"
    )
    checked = 0
    for surface_id, path in sources.items():
        for match in unit_pattern.finditer(path.read_text(encoding="utf-8")):
            unit = match.group(1)
            if not re.search(r"(?im)^\s*entry\b", unit):
                continue
            if not re.search(r"(?im)^\s*save\s*$", unit):
                raise SystemExit(f"{surface_id}: ENTRY unit lacks blanket SAVE")
            checked += 1
    if checked == 0:
        raise SystemExit("corpus SAVE check did not inspect any ENTRY units")


if __name__ == "__main__":
    test_internal_coordinate_order()
    test_save_entry_locals_function_headers()
    test_all_generated_entry_units_have_save()
    print("internal-coordinate order and corpus SAVE checks passed")
