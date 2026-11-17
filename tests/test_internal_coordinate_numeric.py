#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Numeric regression for triatomic internal-coordinate wrapper ordering."""

from __future__ import annotations

import ctypes
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path


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


def _load_library(build_root: Path) -> ctypes.CDLL:
    candidates = sorted(build_root.glob("libpotlibc.so*"))
    if not candidates:
        raise SystemExit("cannot locate libpotlibc.so in build tree")
    return ctypes.CDLL(str(candidates[0]))


def _api_setup(api: ctypes.CDLL) -> None:
    api.potlibc_create.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p)]
    api.potlibc_create.restype = ctypes.c_int
    api.potlibc_destroy.argtypes = [ctypes.c_void_p]
    api.potlibc_destroy.restype = None
    api.potlibc_get_metadata.argtypes = [ctypes.c_void_p, ctypes.POINTER(PotlibcMetadata)]
    api.potlibc_get_metadata.restype = ctypes.c_int
    api.potlibc_evaluate.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    api.potlibc_evaluate.restype = ctypes.c_int
    api.potlibc_energy_directional_derivative.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    api.potlibc_energy_directional_derivative.restype = ctypes.c_int
    api.potlibc_hessian_vector_product.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    api.potlibc_hessian_vector_product.restype = ctypes.c_int


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def _fortran_reference(build_root: Path, surface_id: str, r12: float, r23: float, r13: float) -> tuple[float, float]:
    def fmt_d(value: float) -> str:
        return f"{value:.17E}".replace("E", "D")

    program = f"""\
      program eval
      use dual2_mod
      implicit double precision (a-h,o-z)
      type(dual2) r(3),energy,dedr(3),ec,ew
      integer nsurf,nder,ndum(21)
      type(dual2) easyab,easybc,easyac
      common /zpt1cm/ r,energy,dedr
      common /zpt2cm/ nsurf,nder,ndum
      common /zpt5cm/ easyab,easybc,easyac
      nsurf=0
      nder=1
      call zprepot
      r(1)=dual2({fmt_d(r12)},0d0,0d0,0d0)
      r(2)=dual2({fmt_d(r23)},0d0,0d0,0d0)
      r(3)=dual2({fmt_d(r13)},0d0,0d0,0d0)
      call zpot
      ec=energy
      r(1)=dual2({fmt_d(r12)},0d0,0d0,0d0)
      r(2)=dual2({fmt_d(r13)},0d0,0d0,0d0)
      r(3)=dual2({fmt_d(r23)},0d0,0d0,0d0)
      call zpot
      ew=energy
      write(*,'(a,1x,es24.16)') 'CORRECT',ec%re
      write(*,'(a,1x,es24.16)') 'WRONG',ew%re
      end
"""
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        driver = work / "driver.f"
        binary = work / "driver.exe"
        driver.write_text(program, encoding="ascii")
        surface = build_root / f"{surface_id}-surface.f90"
        utility = build_root / f"{surface_id}-utility.f90"
        dual2 = Path(__file__).resolve().parents[1] / "src" / "dual2_mod.f90"
        compile_run = subprocess.run(
            [
                "gfortran",
                "-ffree-line-length-none",
                str(dual2),
                str(utility),
                str(surface),
                str(driver),
                "-llapack",
                "-lblas",
                "-o",
                str(binary),
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if compile_run.returncode != 0:
            raise SystemExit(
                f"Fortran compile failed for {surface_id}:\n{compile_run.stdout}"
            )
        out = subprocess.run(
            [str(binary)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        ).stdout
    m_correct = re.search(r"CORRECT\s+([+-]?\d+\.\d+E[+-]\d+)", out)
    m_wrong = re.search(r"WRONG\s+([+-]?\d+\.\d+E[+-]\d+)", out)
    if m_correct is None or m_wrong is None:
        raise SystemExit("cannot parse Fortran reference output")
    return float(m_correct.group(1)), float(m_wrong.group(1))


def _check_dual_consistency(
    api: ctypes.CDLL, handle: ctypes.c_void_p, coords: list[float], natoms: int
) -> None:
    n = 3 * natoms
    force = (ctypes.c_double * n)()
    energy = (ctypes.c_double * 1)()
    x = (ctypes.c_double * n)(*coords)
    if api.potlibc_evaluate(handle, 1, x, energy, force) != 0:
        raise SystemExit("potlibc_evaluate failed")
    direction = [((i % 7) - 3) / 5.0 for i in range(n)]
    d = (ctypes.c_double * n)(*direction)
    directional = (ctypes.c_double * 1)()
    if api.potlibc_energy_directional_derivative(handle, 1, x, d, directional) != 0:
        raise SystemExit("directional derivative call failed")
    projected_force = -sum(force[i] * direction[i] for i in range(n))
    if abs(directional[0] - projected_force) > 1.0e-8:
        raise SystemExit("directional derivative does not match projected gradient")

    hvp = (ctypes.c_double * n)()
    if api.potlibc_hessian_vector_product(handle, 1, x, d, hvp) != 0:
        raise SystemExit("hvp call failed")
    step = 1.0e-4
    x_plus = (ctypes.c_double * n)(*[coords[i] + step * direction[i] for i in range(n)])
    x_minus = (ctypes.c_double * n)(*[coords[i] - step * direction[i] for i in range(n)])
    deriv_plus = (ctypes.c_double * 1)()
    deriv_minus = (ctypes.c_double * 1)()
    if api.potlibc_energy_directional_derivative(handle, 1, x_plus, d, deriv_plus) != 0:
        raise SystemExit("displaced + directional derivative failed")
    if api.potlibc_energy_directional_derivative(handle, 1, x_minus, d, deriv_minus) != 0:
        raise SystemExit("displaced - directional derivative failed")
    fd = (deriv_plus[0] - deriv_minus[0]) / (2.0 * step)
    proj = sum(hvp[i] * direction[i] for i in range(n))
    if abs(proj - fd) > 3.0e-5:
        raise SystemExit("HVP projection does not match finite difference")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    build_root = Path(os.environ["POTLIBC_PLUGIN_DIR"]).resolve()
    api = _load_library(build_root)
    _api_setup(api)
    surfaces = [
        (
            "fh25sec",
            # F-H-H nearly collinear transition-state-like geometry
            [(0.0, 0.0, 0.0), (2.6844, 0.0, 0.0), (4.4434, 0.0, 0.0)],
        ),
        (
            "brh2leps",
            # Br-H-H nearly collinear transition-state-like geometry
            [(0.0, 0.0, 0.0), (3.0802, 0.0, 0.0), (4.8202, 0.0, 0.0)],
        ),
        (
            "brhcleleps2",
            # Br-H-Cl bent transition-state-like geometry
            [(0.0, 0.0, 0.0), (3.1, 0.0, 0.0), (1.0, 2.2, 0.0)],
        ),
    ]
    for surface_id, geometry in surfaces:
        handle = ctypes.c_void_p()
        if api.potlibc_create(surface_id.encode(), ctypes.byref(handle)) != 0:
            raise SystemExit(f"potlibc_create failed for {surface_id}")
        metadata = PotlibcMetadata()
        if api.potlibc_get_metadata(handle, ctypes.byref(metadata)) != 0:
            raise SystemExit(f"metadata failed for {surface_id}")
        if metadata.atom_count != 3:
            raise SystemExit(f"{surface_id}: expected 3 atoms")
        xyz = [component for atom in geometry for component in atom]
        x = (ctypes.c_double * 9)(*xyz)
        energy = (ctypes.c_double * 1)()
        force = (ctypes.c_double * 9)()
        if api.potlibc_evaluate(handle, 1, x, energy, force) != 0:
            raise SystemExit(f"evaluate failed for {surface_id}")
        r12 = _distance(geometry[0], geometry[1])
        r23 = _distance(geometry[1], geometry[2])
        r13 = _distance(geometry[0], geometry[2])
        reference, wrong = _fortran_reference(build_root, surface_id, r12, r23, r13)
        if abs(energy[0] - reference) > 1.0e-10:
            raise SystemExit(
                f"{surface_id}: wrapper energy {energy[0]:.16e} "
                f"!= Fortran header-order energy {reference:.16e}"
            )
        if abs(energy[0] - wrong) < 1.0e-6:
            raise SystemExit(
                f"{surface_id}: wrapper still matches legacy AB/AC/BC order"
            )
        _check_dual_consistency(api, handle, xyz, 3)
        api.potlibc_destroy(handle)
    print("numeric internal-coordinate and dual-consistency checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
