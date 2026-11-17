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

import re
import subprocess
import sys
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[1]
generator = root / "tools" / "generate_surface.py"
vendor = root / "upstream" / "potlib2023"

surface_expectations = {
    # Header definitions:
    #   R(1) = R(F-H),   R(2) = R(H-H),   R(3) = R(H-F)
    "fh25sec": ("(1,:)", "(2,:)", "(3,:)"),
    #   R(1) = R(Br-H),  R(2) = R(H-H),   R(3) = R(Br-H)
    "brh2leps": ("(1,:)", "(2,:)", "(3,:)"),
    #   R(1) = R(H-Cl),  R(2) = R(H-H'),  R(3) = R(H'-Cl)
    "clh2g3": ("(1,:)", "(2,:)", "(3,:)"),
}

for surface_id in surface_expectations:
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        subprocess.run(
            [
                sys.executable,
                str(generator),
                str(vendor),
                surface_id,
                str(temp / "source.f90"),
                str(temp / "utility.f90"),
                str(temp / "adapter.f90"),
            ],
            check=True,
        )
        adapter = (temp / "adapter.f90").read_text(encoding="utf-8")
        assignments = re.findall(
            r"r\((\d)\)=distance3\(coords\((\d),:\),coords\((\d),:\)\)",
            adapter,
        )
        if assignments[:3] != [("1", "1", "2"), ("2", "2", "3"), ("3", "1", "3")]:
            raise SystemExit(
                f"{surface_id}: expected R order (AB,BC,AC) from header, "
                f"got {assignments[:3]!r}"
            )

print("internal-coordinate order checks passed for FH2, BrH2, and ClH2")
