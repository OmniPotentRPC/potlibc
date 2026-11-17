#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Build the documentation with every warning promoted to an error."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "docs" / "source"


def run(*arguments: str) -> None:
    subprocess.run(arguments, cwd=ROOT, check=True)


def main() -> int:
    shutil.rmtree(SOURCE / "_doxygen", ignore_errors=True)
    shutil.rmtree(ROOT / "docs" / "_build", ignore_errors=True)
    run(
        sys.executable,
        "docs/generate_surface_catalog.py",
        "data/surfaces.json",
        "docs/source/surfaces.rst",
    )
    run("doxygen", "docs/Doxyfile")
    run(
        sys.executable,
        "-m",
        "sphinx",
        "-W",
        "--keep-going",
        "-d",
        "docs/_build/doctrees",
        "docs/source",
        "docs/_build/html",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
