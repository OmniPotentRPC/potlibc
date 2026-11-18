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
import sys
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


def _fortran_reference(
    build_root: Path,
    surface_id: str,
    distances: tuple[float, float, float],
    geometry: list[tuple[float, float, float]],
) -> float:
    def fmt_d(value: float) -> str:
        return f"{value:.17E}".replace("E", "D")

    if surface_id in {"hf2jotii2001", "hi2kka2001"}:
        assignments = "\n".join(
            f"      cart({atom},{axis})=dual2({fmt_d(value)},0d0,0d0,0d0)"
            for atom, coordinates in enumerate(geometry, 1)
            for axis, value in enumerate(coordinates, 1)
        )
        program = f"""\
      program eval
      use dual2_mod
      implicit none
      type(dual2) r(75),energy,dedr(75),cart(25,3),anuzero
      integer nullbl(25),nflag(20),nasurf(6,6),nder
      common /zpt1cm/ r,energy,dedr
      common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
      call zprepot
{assignments}
      call zpot
      write(*,'(a,1x,es24.16)') 'ENERGY',energy%re
      end
"""
    elif surface_id in {"nafh2d", "nafh2v"}:
        program = f"""\
      program eval
      use dual2_mod
      implicit none
      type(dual2) r(1,3),energy(1),dedr(3,1)
      integer g_myrank,g_nprocs
      common /zcom_para/ g_myrank,g_nprocs
      g_myrank=1
      g_nprocs=1
      call zprepot
      r(1,1)=dual2({fmt_d(distances[0])},0d0,0d0,0d0)
      r(1,2)=dual2({fmt_d(distances[1])},0d0,0d0,0d0)
      r(1,3)=dual2({fmt_d(distances[2])},0d0,0d0,0d0)
      call zpot(r,energy,dedr,1,1)
      write(*,'(a,1x,es24.16)') 'ENERGY',energy(1)%re
      end
"""
    else:
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
      r(1)=dual2({fmt_d(distances[0])},0d0,0d0,0d0)
      r(2)=dual2({fmt_d(distances[1])},0d0,0d0,0d0)
      r(3)=dual2({fmt_d(distances[2])},0d0,0d0,0d0)
      call zpot
      ec=energy
      write(*,'(a,1x,es24.16)') 'ENERGY',ec%re
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
    match = re.search(r"ENERGY\s+([+-]?\d+\.\d+E[+-]\d+)", out)
    if match is None:
        raise SystemExit("cannot parse Fortran reference output")
    return float(match.group(1))


def _evaluate_energy(
    api: ctypes.CDLL, surface_id: str, coords: list[float]
) -> float:
    handle = ctypes.c_void_p()
    if api.potlibc_create(surface_id.encode(), ctypes.byref(handle)) != 0:
        raise SystemExit(f"potlibc_create failed for {surface_id}")
    try:
        xyz = (ctypes.c_double * len(coords))(*coords)
        energy = (ctypes.c_double * 1)()
        force = (ctypes.c_double * len(coords))()
        if api.potlibc_evaluate(handle, 1, xyz, energy, force) != 0:
            raise SystemExit(f"potlibc_evaluate failed for {surface_id}")
        return energy[0]
    finally:
        api.potlibc_destroy(handle)


def _check_repeated_atom_symmetry(
    api: ctypes.CDLL,
) -> None:
    geometry = [
        0.0, 0.0, 0.0,
        2.3, 0.2, -0.1,
        0.4, 1.9, 0.3,
    ]
    failures: list[str] = []
    cases = {
        "hi2kka": (1, 2),
        "hf2jotii": (0, 1),
        "hf2jotii2001": (1, 2),
        "hi2kka2001": (0, 2),
    }
    for surface_id, (left, right) in cases.items():
        original = _evaluate_energy(api, surface_id, geometry)
        swapped = geometry.copy()
        swapped[3 * left:3 * left + 3], swapped[3 * right:3 * right + 3] = (
            swapped[3 * right:3 * right + 3],
            swapped[3 * left:3 * left + 3],
        )
        permuted = _evaluate_energy(api, surface_id, swapped)
        scale = max(1.0, abs(original), abs(permuted))
        tolerance = 64.0 * sys.float_info.epsilon * scale
        if (
            not math.isfinite(original)
            or not math.isfinite(permuted)
            or abs(original - permuted) > tolerance
        ):
            failures.append(
                f"{surface_id} atoms {left + 1}/{right + 1}: "
                f"{original:.16e} != {permuted:.16e} "
                f"(tol {tolerance:.3e})"
            )
    if failures:
        raise SystemExit(
            "repeated-atom symmetry failures:\n" + "\n".join(failures)
        )


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


def _check_energy_fd(
    api: ctypes.CDLL,
    handle: ctypes.c_void_p,
    coords: list[float],
    direction: list[float],
) -> tuple[float, float, float]:
    n = len(coords)
    x = (ctypes.c_double * n)(*coords)
    d = (ctypes.c_double * n)(*direction)
    energy = (ctypes.c_double * 1)()
    force = (ctypes.c_double * n)()
    directional = (ctypes.c_double * 1)()
    if api.potlibc_evaluate(handle, 1, x, energy, force) != 0:
        raise SystemExit("potlibc_evaluate failed")
    if api.potlibc_energy_directional_derivative(handle, 1, x, d, directional) != 0:
        raise SystemExit("directional derivative call failed")
    grad_proj = -sum(force[i] * direction[i] for i in range(n))
    if abs(directional[0] - grad_proj) > 1.0e-10:
        raise SystemExit(
            f"directional derivative {directional[0]:.16e} != "
            f"gradient projection {grad_proj:.16e}"
        )
    for step, tol in ((1.0e-4, 5.0e-7), (1.0e-5, 5.0e-8), (1.0e-6, 5.0e-7)):
        xp = (ctypes.c_double * n)(*[coords[i] + step * direction[i] for i in range(n)])
        xm = (ctypes.c_double * n)(*[coords[i] - step * direction[i] for i in range(n)])
        ep = (ctypes.c_double * 1)()
        em = (ctypes.c_double * 1)()
        fp = (ctypes.c_double * n)()
        fm = (ctypes.c_double * n)()
        if api.potlibc_evaluate(handle, 1, xp, ep, fp) != 0:
            raise SystemExit("displaced + evaluate failed")
        if api.potlibc_evaluate(handle, 1, xm, em, fm) != 0:
            raise SystemExit("displaced - evaluate failed")
        fd = (ep[0] - em[0]) / (2.0 * step)
        if abs(fd - directional[0]) > tol or abs(fd - grad_proj) > tol:
            raise SystemExit(
                f"step {step}: central FD {fd:.16e} disagrees with "
                f"directional {directional[0]:.16e} / grad_proj {grad_proj:.16e}"
            )
    return energy[0], directional[0], grad_proj


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    build_root = Path(
        os.environ.get("POTLIBC_PLUGIN_DIR", str(root / "build"))
    ).resolve()
    os.environ.setdefault("POTLIBC_PLUGIN_DIR", str(build_root))
    api = _load_library(build_root)
    _api_setup(api)
    surfaces = [
        (
            "fh25sec",
            # F-H-H nearly collinear transition-state-like geometry
            [(0.0, 0.0, 0.0), (2.6844, 0.0, 0.0), (4.4434, 0.0, 0.0)],
            ("r12", "r23", "r13"),
        ),
        (
            "brh2leps",
            # Br-H-H nearly collinear transition-state-like geometry
            [(0.0, 0.0, 0.0), (3.0802, 0.0, 0.0), (4.8202, 0.0, 0.0)],
            ("r12", "r23", "r13"),
        ),
        (
            "brhcleleps2",
            # Br-H-Cl bent transition-state-like geometry
            [(0.0, 0.0, 0.0), (3.1, 0.0, 0.0), (1.0, 2.2, 0.0)],
            ("r12", "r23", "r13"),
        ),
        (
            "clh2g3",
            # Cl-H-H bent encounter geometry
            [(0.0, 0.0, 0.0), (2.4, 0.0, 0.0), (0.2, 1.5, 0.1)],
            ("r12", "r23", "r13"),
        ),
        (
            "clh2alab",
            # Cl-H-H bent encounter geometry
            [(0.0, 0.0, 0.0), (2.4, 0.0, 0.0), (0.2, 1.8, 0.1)],
            ("r12", "r23", "r13"),
        ),
        (
            "hi2kka",
            [(0.0, 0.0, 0.0), (3.0, 0.1, 0.0), (4.8, 1.2, 0.2)],
            ("r12", "r13", "r23"),
        ),
        (
            "hf2jotii",
            [(0.0, 0.0, 0.0), (2.7, 0.2, 0.0), (1.0, 1.6, 0.1)],
            ("r13", "r12", "r23"),
        ),
        (
            "nafh2d",
            [(0.0, 0.0, 0.0), (3.0, 0.2, 0.0), (4.2, 1.4, 0.1)],
            ("r12", "r23", "r13"),
        ),
        (
            "nafh2v",
            [(0.0, 0.0, 0.0), (3.0, 0.2, 0.0), (4.2, 1.4, 0.1)],
            ("r12", "r23", "r13"),
        ),
        (
            "hf2jotii2001",
            [(0.0, 0.0, 0.0), (1.7, 0.2, 0.0), (3.5, 1.2, 0.1)],
            ("r12", "r23", "r13"),
        ),
        (
            "hi2kka2001",
            [(0.0, 0.0, 0.0), (3.0, 0.1, 0.0), (4.8, 1.2, 0.2)],
            ("r12", "r23", "r13"),
        ),
    ]
    source_atom_numbers = {
        "nafh2d": [11, 1, 9],
        "nafh2v": [11, 1, 9],
        "hf2jotii2001": [1, 9, 9],
        "hi2kka2001": [53, 1, 53],
    }
    coordinate_failures: list[str] = []
    for surface_id, geometry, coordinate_order in surfaces:
        handle = ctypes.c_void_p()
        if api.potlibc_create(surface_id.encode(), ctypes.byref(handle)) != 0:
            raise SystemExit(f"potlibc_create failed for {surface_id}")
        metadata = PotlibcMetadata()
        if api.potlibc_get_metadata(handle, ctypes.byref(metadata)) != 0:
            raise SystemExit(f"metadata failed for {surface_id}")
        if metadata.atom_count != 3:
            raise SystemExit(f"{surface_id}: expected 3 atoms")
        wrapper_geometry = geometry
        if surface_id in source_atom_numbers:
            source_numbers = source_atom_numbers[surface_id]
            unused = list(range(3))
            wrapper_geometry = []
            for index in range(3):
                number = metadata.atomic_numbers[index]
                match = next(
                    (
                        source_index
                        for source_index in unused
                        if source_numbers[source_index] == number
                    ),
                    None,
                )
                if match is None:
                    raise SystemExit(
                        f"{surface_id}: published atom {number} is not in source order"
                    )
                unused.remove(match)
                wrapper_geometry.append(geometry[match])
        xyz = [component for atom in wrapper_geometry for component in atom]
        x = (ctypes.c_double * 9)(*xyz)
        energy = (ctypes.c_double * 1)()
        force = (ctypes.c_double * 9)()
        if api.potlibc_evaluate(handle, 1, x, energy, force) != 0:
            raise SystemExit(f"evaluate failed for {surface_id}")
        r12 = _distance(geometry[0], geometry[1])
        r23 = _distance(geometry[1], geometry[2])
        r13 = _distance(geometry[0], geometry[2])
        distances = {"r12": r12, "r23": r23, "r13": r13}
        reference = _fortran_reference(
            build_root,
            surface_id,
            tuple(distances[name] for name in coordinate_order),
            geometry,
        )
        if abs(energy[0] - reference) > 1.0e-10:
            coordinate_failures.append(
                f"{surface_id}: wrapper energy {energy[0]:.16e} "
                f"!= Fortran header-order energy {reference:.16e}"
            )
        _check_dual_consistency(api, handle, xyz, 3)
        api.potlibc_destroy(handle)
    if coordinate_failures:
        raise SystemExit(
            "coordinate reference failures:\n" + "\n".join(coordinate_failures)
        )

    clh2_geometries = [
        [
            0.0, 0.0, 0.0,
            1.0 + 1.0 / 31.0, 1.0 / 47.0, 1.0 / 59.0,
            2.0 / 31.0, 1.0 + 4.0 / 47.0, 8.0 / 59.0,
        ],
        [0.0, 0.0, 0.0, 2.4, 0.0, 0.0, 0.2, 1.5, 0.1],
        [0.1, -0.2, 0.05, 2.7, 0.3, -0.1, 2.9, 1.7, 0.2],
    ]
    direction = [(i + 1) / 9.0 for i in range(9)]
    h_g3 = ctypes.c_void_p()
    h_g3v = ctypes.c_void_p()
    if api.potlibc_create(b"clh2g3", ctypes.byref(h_g3)) != 0:
        raise SystemExit("potlibc_create failed for clh2g3")
    if api.potlibc_create(b"clh2g3v", ctypes.byref(h_g3v)) != 0:
        raise SystemExit("potlibc_create failed for clh2g3v")
    for coords in clh2_geometries:
        e_g3, d_g3, _ = _check_energy_fd(api, h_g3, coords, direction)
        e_g3v, d_g3v, _ = _check_energy_fd(api, h_g3v, coords, direction)
        if abs(e_g3 - e_g3v) > 5.0e-6 or abs(d_g3 - d_g3v) > 5.0e-6:
            raise SystemExit(
                f"clh2g3 and clh2g3v disagree: ({e_g3}, {d_g3}) vs ({e_g3v}, {d_g3v})"
            )
    api.potlibc_destroy(h_g3)
    api.potlibc_destroy(h_g3v)

    _check_repeated_atom_symmetry(api)

    print("numeric internal-coordinate and dual-consistency checks passed")
    return 0


def test_internal_coordinate_numeric() -> None:
    assert main() == 0


if __name__ == "__main__":
    raise SystemExit(main())
