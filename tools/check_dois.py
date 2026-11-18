#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Verify every DOI cited by potlibc resolves through the handle API.

This checks the Handle System API (the same one ``doi.org`` uses to
redirect) directly, so it is not affected by publisher sites that return
4xx responses to scripted HTTP clients. See ``docs/source/accuracy.rst``
and ``sphinx`` ``linkcheck`` for the companion link check, which covers
plain HTTP(S) URLs instead of DOIs.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable

HANDLE_API = "https://doi.org/api/handles/{doi}"

# Repository DOIs that are not recorded in data/surfaces.json because they
# describe the vendor deposit itself rather than one surface.
REPOSITORY_DOIS = ["10.17632/chxxhnnt4p.1"]


def surface_dois(catalog_path: Path) -> list[str]:
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    dois: set[str] = set()
    for surface in catalog["surfaces"]:
        dois.update(surface["documentation"]["dois"])
    return sorted(dois)


def check_doi(doi: str) -> tuple[bool, str]:
    url = HANDLE_API.format(doi=doi)
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            payload = json.loads(response.read())
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as error:
        return False, f"request failed: {error}"
    response_code = payload.get("responseCode")
    if response_code != 1:
        return False, f"handle lookup failed (responseCode={response_code})"
    return True, "resolved"


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "catalog",
        type=Path,
        help="Path to data/surfaces.json",
    )
    args = parser.parse_args(argv)

    dois = sorted(set(surface_dois(args.catalog)) | set(REPOSITORY_DOIS))
    failures = []
    for doi in dois:
        ok, detail = check_doi(doi)
        status = "OK" if ok else "FAIL"
        print(f"{status} {doi}: {detail}")
        if not ok:
            failures.append(doi)

    print(f"Checked {len(dois)} DOIs, {len(failures)} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
