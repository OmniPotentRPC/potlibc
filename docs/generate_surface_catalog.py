#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Generate the exhaustive Sphinx surface catalog from surfaces.json."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def rst(value: str) -> str:
    """Escape deposited text for a single reStructuredText paragraph."""
    return (
        " ".join(value.split())
        .replace("\\", "\\\\")
        .replace("*", "\\*")
        .replace("`", "\\`")
        .replace("|", "\\|")
    )


def atom_count(value: dict[str, Any]) -> str:
    if value["kind"] == "fixed":
        return str(value["value"])
    if value["kind"] == "variable":
        return f"variable ({value['symbol']})"
    return "not determined"


def missing(label: str) -> str:
    return f"Not recorded in the deposited {label}."


def render(catalog: dict[str, Any]) -> str:
    lines = [
        ".. SPDX-FileCopyrightText: 2026 Rohit Goswami",
        ".. SPDX-License-Identifier: BSD-3-Clause",
        "",
        "Surface catalog",
        "===============",
        "",
        "This page is generated from ``data/surfaces.json``. Descriptive text is",
        "copied from each deposited source header; an explicit ``Not recorded``",
        "means that the vendor header contains no corresponding field.",
        "",
        f"The catalog contains **{len(catalog['surfaces'])}** callable surfaces.",
        "",
    ]
    for surface in catalog["surfaces"]:
        documentation = surface["documentation"]
        identifier = surface["id"]
        title = f"{identifier} - {surface['name']}"
        source_url = (
            "https://github.com/OmniPotentRPC/potlibc/blob/main/"
            f"upstream/potlib2023/{surface['path']}"
        )
        lines.extend(
            [
                f".. _surface-{identifier}:",
                "",
                title,
                "-" * len(title),
                "",
                f":Formula / system: {rst(documentation['formula']) or missing('System field')}",
                f":Atoms: {atom_count(surface['atom_count'])}",
                f":Interface: ``{surface['interface_family']}``",
                f":Coordinates: {rst(documentation['coordinates']) or missing('Coordinates field')}",
                f":Original reference: {rst(documentation['reference']) or missing('Reference field')}",
            ]
        )
        if documentation["dois"]:
            links = ", ".join(
                f"`{rst(doi)} <https://doi.org/{doi}>`_"
                for doi in documentation["dois"]
            )
            lines.append(f":DOI: {links}")
        else:
            lines.append(":DOI: Not recorded in the deposited source header.")
        lines.extend(
            [
                f":Source: `{surface['path']} <{source_url}>`_",
                f":SHA-256: ``{surface['sha256']}``",
                "",
            ]
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render(catalog), encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
