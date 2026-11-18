#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Check deposited Angstrom and kcal/mol or eV adapters honor Bohr and Hartree."""

from __future__ import annotations

import ctypes
import math
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import generate_surface as generator  # noqa: E402


BOHR_TO_ANGSTROM = 0.529177210903
KCAL_PER_HARTREE = 627.5094740631
BOHR_TOKEN = "0.529177210903d0"
KCAL_TOKEN = "627.5094740631d0"
EV_TOKEN = "27.211386245988d0"


class PotlibcMetadata(ctypes.Structure):
    _fields_ = [
        ("surface_id", ctypes.c_char_p),
        ("atom_count", ctypes.c_size_t),
        ("atomic_numbers", ctypes.POINTER(ctypes.c_int)),
        ("reference_geometry", ctypes.POINTER(ctypes.c_double)),
        ("atom_ordering", ctypes.c_char_p),
        ("length_unit", ctypes.c_char_p),
        ("energy_unit", ctypes.c_char_p),
        ("force_unit", ctypes.c_char_p),
        ("capabilities", ctypes.c_uint32),
    ]


def _bridges(vendor: Path) -> dict[str, str]:
    catalog = generator.inventory(vendor)
    by_id = {item["id"]: item for item in catalog["surfaces"]}
    wanted = [
        "n3pes",
        "o3pes",
        "o4pes",
        "c2h7",
        "c2h7pot5",
        "no2-2ap-pip-nn",
        "no2-4ap-pip-nn",
        "no2-6ap-pip-nn",
        "fh2o-pipnn",
        "soc-fh2o-nn",
        "ho3-pes",
        "h3o2-pipnn",
        "pes-h3cl",
        "ch3nh2",
        "pes-shell",
        "phsch3-ia",
        "c2h6cl-rev2019",
        "pes-n4-mb-pip-nn",
    ]
    return {
        surface_id: generator._energy_bridge(
            generator._plan(by_id[surface_id], vendor)
        )
        for surface_id in wanted
    }


def _require(text: str, token: str, label: str, failures: list[str]) -> None:
    if token not in text:
        failures.append(f"{label}: missing {token}")


def _forbid(text: str, token: str, label: str, failures: list[str]) -> None:
    if token in text:
        failures.append(f"{label}: unexpected {token}")


def _check_bridges(bridges: dict[str, str]) -> None:
    failures: list[str] = []
    for surface_id in ("n3pes", "o3pes"):
        text = bridges[surface_id]
        _require(text, BOHR_TOKEN, surface_id, failures)
        _require(text, f"value=z{surface_id}(r12,r13,r23)/{KCAL_TOKEN}", surface_id, failures)
    text = bridges["o4pes"]
    _require(text, f"r=r*{BOHR_TOKEN}", "o4pes", failures)
    _require(text, f"value=zo4pes(r)/{KCAL_TOKEN}", "o4pes", failures)
    for surface_id in ("c2h7", "c2h7pot5"):
        text = bridges[surface_id]
        _require(text, f"xyz=coords*{BOHR_TOKEN}", surface_id, failures)
        _require(text, "call zc2h7pot(xyz,value)", surface_id, failures)
        _require(text, f"value=value/{KCAL_TOKEN}", surface_id, failures)
        _forbid(text, "reshape(transpose", surface_id, failures)
    for surface_id in (
        "no2-2ap-pip-nn",
        "no2-4ap-pip-nn",
        "no2-6ap-pip-nn",
        "fh2o-pipnn",
        "ho3-pes",
        "h3o2-pipnn",
        "pes-h3cl",
    ):
        text = bridges[surface_id]
        _require(text, BOHR_TOKEN, surface_id, failures)
        _require(text, f"value=value/{EV_TOKEN}", surface_id, failures)
        _forbid(text, KCAL_TOKEN, surface_id, failures)
    text = bridges["soc-fh2o-nn"]
    _require(text, BOHR_TOKEN, "soc-fh2o-nn", failures)
    _require(text, f"value=vg/{KCAL_TOKEN}", "soc-fh2o-nn", failures)
    _forbid(text, EV_TOKEN, "soc-fh2o-nn", failures)
    text = bridges["ch3nh2"]
    _require(text, f"xx=transpose(coords)*{BOHR_TOKEN}", "ch3nh2", failures)
    _forbid(text, KCAL_TOKEN, "ch3nh2", failures)
    _forbid(text, EV_TOKEN, "ch3nh2", failures)
    text = bridges["pes-shell"]
    _require(text, f"matrix=transpose(coords)*{BOHR_TOKEN}", "pes-shell", failures)
    _require(text, "value=zf(matrix)", "pes-shell", failures)
    _forbid(text, KCAL_TOKEN, "pes-shell", failures)
    _forbid(text, EV_TOKEN, "pes-shell", failures)
    text = bridges["phsch3-ia"]
    _require(text, "xx=transpose(coords)", "phsch3-ia", failures)
    _forbid(text, BOHR_TOKEN, "phsch3-ia", failures)
    text = bridges["c2h6cl-rev2019"]
    _require(text, f"q(k)=coords(atom,axis)*{BOHR_TOKEN}", "c2h6cl-rev2019", failures)
    _require(text, "value=pengygs", "c2h6cl-rev2019", failures)
    _forbid(text, KCAL_TOKEN, "c2h6cl-rev2019", failures)
    _forbid(text, EV_TOKEN, "c2h6cl-rev2019", failures)
    text = bridges["pes-n4-mb-pip-nn"]
    _require(text, "call zn4pipnn", "pes-n4-mb-pip-nn", failures)
    _forbid(text, BOHR_TOKEN, "pes-n4-mb-pip-nn", failures)
    if failures:
        raise SystemExit("adapter unit mismatches:\n" + "\n".join(failures))


def _n3_reference_hartree(vendor: Path) -> float:
    terms: list[tuple[float, float, float, float]] = []
    coeff_path = vendor / "src" / "done" / "N3" / "coeff_n3_o3.txt"
    for line in coeff_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split()
        terms.append(
            (float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]))
        )
    if len(terms) != 52:
        raise SystemExit(f"expected 52 N3 coefficients, found {len(terms)}")
    re = 1.098
    ae = 0.9
    distances = (1.098, 2.0, 2.0)
    morse = tuple(math.exp(-(distance - re) / ae) for distance in distances)
    energy = 0.0
    for i1, i2, i3, coeff in terms:
        x, y, z = morse
        energy += coeff * (
            x**i1 * y**i2 * z**i3
            + x**i1 * y**i3 * z**i2
            + x**i2 * y**i3 * z**i1
            + x**i2 * y**i1 * z**i3
            + x**i3 * y**i1 * z**i2
            + x**i3 * y**i2 * z**i1
        )
    return energy / KCAL_PER_HARTREE


def _n3_geometry() -> list[float]:
    r12, r13, r23 = (1.098 / BOHR_TO_ANGSTROM, 2.0 / BOHR_TO_ANGSTROM, 2.0 / BOHR_TO_ANGSTROM)
    cx = (r13 * r13 + r12 * r12 - r23 * r23) / (2.0 * r12)
    cy = math.sqrt(r13 * r13 - cx * cx)
    return [0.0, 0.0, 0.0, r12, 0.0, 0.0, cx, cy, 0.0]


def _check_n3_energy(vendor: Path, build_root: Path) -> None:
    candidates = sorted(build_root.glob("libpotlibc.so*"))
    if not candidates:
        raise SystemExit(f"cannot locate libpotlibc.so in {build_root}")
    api = ctypes.CDLL(str(candidates[0]))
    api.potlibc_create.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p)]
    api.potlibc_create.restype = ctypes.c_int
    api.potlibc_destroy.argtypes = [ctypes.c_void_p]
    api.potlibc_destroy.restype = None
    api.potlibc_get_metadata.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(PotlibcMetadata),
    ]
    api.potlibc_get_metadata.restype = ctypes.c_int
    api.potlibc_evaluate.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    api.potlibc_evaluate.restype = ctypes.c_int
    handle = ctypes.c_void_p()
    if api.potlibc_create(b"n3pes", ctypes.byref(handle)) != 0:
        raise SystemExit("potlibc_create failed for n3pes")
    try:
        metadata = PotlibcMetadata()
        if api.potlibc_get_metadata(handle, ctypes.byref(metadata)) != 0:
            raise SystemExit("potlibc_get_metadata failed for n3pes")
        if metadata.energy_unit != b"hartree" or metadata.length_unit != b"bohr":
            raise SystemExit(
                "n3pes metadata units are "
                f"{metadata.length_unit!r} and {metadata.energy_unit!r}"
            )
        coords = _n3_geometry()
        xyz = (ctypes.c_double * len(coords))(*coords)
        energy = (ctypes.c_double * 1)()
        force = (ctypes.c_double * len(coords))()
        if api.potlibc_evaluate(handle, 1, xyz, energy, force) != 0:
            raise SystemExit("potlibc_evaluate failed for n3pes")
        expected = _n3_reference_hartree(vendor)
        if abs(energy[0] - expected) > 1.0e-8:
            raise SystemExit(
                f"n3pes energy {energy[0]:.16e} != reference {expected:.16e}"
            )
    finally:
        api.potlibc_destroy(handle)


def main() -> int:
    vendor = ROOT / "upstream" / "potlib2023"
    _check_bridges(_bridges(vendor))
    build_root = Path(os.environ.get("POTLIBC_PLUGIN_DIR", str(ROOT / "build"))).resolve()
    os.environ.setdefault("POTLIBC_PLUGIN_DIR", str(build_root))
    _check_n3_energy(vendor, build_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
