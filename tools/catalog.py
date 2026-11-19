#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Inventory POTLIB sources and generate adapters without changing vendor files."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True

FORTRAN_SUFFIXES = {".f", ".F", ".f90", ".F90"}
HEADER = re.compile(r"(?im)^[c*!]\s*([^:\n]+):\s*(.*?)\s*$")
ACTIVE_CALL = r"(?im)^(?![c*!])([ \t]*call\s+{}\s*(?:\([^)]*\))?[^\n]*)$"


def fields(text: str) -> dict[str, str]:
    return {key.strip().lower(): value.strip() for key, value in HEADER.findall(text[:12000])}


def surface_id(path: Path) -> str:
    return re.sub(r"[^a-z0-9]+", "-", path.stem.lower()).strip("-")


def classify(path: Path, root: Path) -> dict[str, object]:
    relative = path.relative_to(root).as_posix()
    text = path.read_text(errors="replace")
    info = fields(text)
    interface = info.get("interface", "")
    entry: dict[str, object] = {
        "path": relative,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "system": info.get("system", ""),
        "name": info.get("common name", path.stem),
        "interface": interface,
        "wrapped": False,
        "reason": "",
    }
    has_pot = re.search(r"(?im)^\s*subroutine\s+pot\b", text) is not None
    has_prepot = re.search(r"(?im)^\s*subroutine\s+prepot\b", text) is not None
    if has_pot and has_prepot:
        atoms_match = (
            re.search(r"(?im)number\s+of\s+bodies\s*:\s*(\d+)", text[:12000])
            or re.search(r"(?im)\bdata\s+natoms\s*/\s*(\d+)\s*/", text)
        )
        derivatives_match = re.search(
            r"(?im)number\s+of\s+derivatives\s*:\s*(\d+)", text[:12000]
        )
        identifier = surface_id(path)
        entry.update(
            {
                "id": identifier,
                "atoms": int(atoms_match.group(1)) if atoms_match else None,
                "derivatives": (
                    int(derivatives_match.group(1))
                    if derivatives_match
                    else None
                ),
                "wrapped": True,
                "reason": "mechanically adapted canonical POT/PREPOT surface",
                "fixed_unit_init": (
                    len(re.findall(ACTIVE_CALL.format("potinfo"), text)) == 1
                    and len(re.findall(ACTIVE_CALL.format("ancvrt"), text)) == 1
                ),
                "hvp": (
                    identifier == "ch52008"
                    and atoms_match is not None
                    and len(re.findall(ACTIVE_CALL.format("potinfo"), text)) == 1
                    and len(re.findall(ACTIVE_CALL.format("ancvrt"), text)) == 1
                ),
            }
        )
    elif text.startswith("This file contains a bug"):
        entry["reason"] = "vendor marks this revision as buggy and replaced"
    elif re.search(r"(?im)^\s*(program|module)\b", text):
        entry["reason"] = "auxiliary/test/module source, not a standalone surface entry"
    elif info.get("data file") or re.search(
        r"(?im)^\s*[c*!]\s*(?:data|parameter)\s+files?\s*:", text[:12000]
    ):
        entry["reason"] = "requires external parameter files and lacks canonical POT/PREPOT entries"
    elif has_pot:
        entry["reason"] = "POT entry exists without canonical PREPOT initialization"
    elif has_prepot:
        entry["reason"] = "PREPOT entry exists without canonical POT evaluation"
    elif interface:
        entry["reason"] = f"interface family {interface!r} lacks canonical POT/PREPOT entries"
    else:
        entry["reason"] = "no canonical POT/PREPOT surface entries"
    return entry


def inventory(vendor: Path) -> list[dict[str, object]]:
    src = vendor / "src"
    paths = sorted(
        path for path in src.rglob("*") if path.is_file() and path.suffix in FORTRAN_SUFFIXES
    )
    entries = [classify(path, vendor) for path in paths]
    seen: set[str] = set()
    for entry in entries:
        if entry["wrapped"]:
            identifier = str(entry["id"])
            if identifier in seen:
                entry["wrapped"] = False
                entry["reason"] = f"generated surface id {identifier!r} collides with another entry"
                entry.pop("id", None)
            else:
                seen.add(identifier)
    return entries


def nest_variadic_legacy_intrinsics(text: str) -> str:
    """Nest variadic MIN/MAX calls so derived-type overloads stay binary."""
    matches = list(re.finditer(r"(?i)\b(dmin1|dmax1|min|max)\s*\(", text))
    for match in reversed(matches):
        depth = 1
        quote: str | None = None
        arguments: list[str] = []
        start = match.end()
        position = start
        while position < len(text) and depth:
            character = text[position]
            if quote:
                if character == quote:
                    quote = None
            elif character in "'\"":
                quote = character
            elif character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    arguments.append(text[start:position])
                    break
            elif character == "," and depth == 1:
                arguments.append(text[start:position])
                start = position + 1
            position += 1
        if depth or len(arguments) <= 2:
            continue
        name = match.group(1)
        nested = f"{name}({arguments[0]},{arguments[1]})"
        for argument in arguments[2:]:
            nested = f"{name}({nested},{argument})"
        text = text[:match.start()] + nested + text[position + 1:]
    return text


def normalize(text: str) -> str:
    text = nest_variadic_legacy_intrinsics(text)
    text = re.sub(r"(?m)[ \t]+\*[A-Z]{2}\d{6}\s*$", "", text)
    text = re.sub(
        r"(?im)^[ \t]*deep\s+neural\s+network[ \t]*$",
        "! deep neural network",
        text,
    )
    text = re.sub(
        r"(?i)\bisnan\s*\(\s*tmp\s*\(\s*i\s*,\s*j\s*\)\s*\)",
        "tmp(i,j).ne.tmp(i,j)",
        text,
    )
    text = re.sub(
        r"(?i)\bmaxval\s*\(\s*r\s*,\s*nbonds\s*\)",
        "maxval(r)",
        text,
    )
    if re.search(r"(?im)^\s*20\s+COS\s*=", text):
        text = re.sub(r"(?i)\bCOS\b", "COSTERM", text)
    if re.search(r"(?im)^\s*subroutine\s+patchfunc\b", text):
        text = re.sub(
            r"(?ims)^(\s*subroutine\s+patchfunc\b.*?^\s*end\s+subroutine\s+patchfunc\s*$)",
            lambda match: re.sub(r"(?i)\ba\b", "angle", match.group(1)),
            text,
        )
    if (
        re.search(r"(?i)\bcall\s+second\s*\(", text)
        and not re.search(r"(?im)^\s*subroutine\s+second\s*\(", text)
    ):
        text += (
            "\n      SUBROUTINE SECOND(TIME)\n"
            "      DOUBLE PRECISION TIME\n"
            "      TIME=0.0D0\n"
            "      RETURN\n"
            "      END\n"
        )
    text = text.replace(
        "\nand the potential coupling.\n",
        "\nC and the potential coupling.\n",
    )
    if text.startswith("This file contains a bug"):
        lines = text.splitlines(keepends=True)
        lines[0] = "C " + lines[0]
        lines[1] = "C " + lines[1]
        text = "".join(lines)
    if text.startswith("THIS SUBROUTINE IS IDENTICAL"):
        lines = text.splitlines(keepends=True)
        lines[0] = "C " + lines[0]
        lines[1] = "C " + lines[1]
        text = "".join(lines)
    text = re.sub(r"(?im)^[ \t]*Program\s+pesa2\.f\s*$", "! package banner", text)
    text = re.sub(
        r"(?im)^(\s*SUBROUTINE\s+PRPOT2\s*\n\s*IMPLICIT[^\n]*\n)"
        r"(\s*PARAMETER\s*\(\s*HRTREE\s*=)",
        r"\1      REAL*8 HRTREE\n\2",
        text,
    )
    if len(re.findall(r"(?im)^\s*subroutine\s+prepot\b", text)) == 2:
        text = re.sub(
            r"(?im)^(\s*)subroutine\s+prepot\b[^\n]*",
            r"C potlibc removed duplicate PREPOT declaration",
            text,
            count=1,
        )
    if (
        re.search(r"(?i)read\s*\([^)]*\)\s*a\(i\),da\(i\)", text)
        and not re.search(r"(?im)^\s*dimension\s+a\s*\(", text)
    ):
        text = re.sub(
            r"(?im)^(\s*parameter\s*\(\s*ma\s*=\s*135\s*\)[^\n]*)$",
            r"\1\n      dimension a(ma),da(ma)\n      common/par/ a",
            text,
            count=1,
        )
    text = re.sub(
        r"(?im)^(\s*double\s+precision\s+)rhf0did,rhf0disp,"
        r"(hf_did,hf_disp,e_LR_HF\s*)$",
        r"\1\2",
        text,
    )
    text = re.sub(
        r"(?im)^(\s*&\s*coschi,)cevau,(r2LRcut_a,r2LRcut_r\s*)$",
        r"\1\2",
        text,
    )
    if "g_pot11(g_p_, r(i,1)" in text and not re.search(
        r"(?i)dimension\s+v11\s*\(", text
    ):
        text = text.replace(
            "dimension de1(3,1),de2(3,1),de3(3,1),ddis(3)",
            "dimension de1(3,1),de2(3,1),de3(3,1),ddis(3)\n"
            "        dimension v11(1),v12(1),v22(1)",
            1,
        )
        for name in ("v11", "v12", "v22"):
            text = re.sub(
                rf"\b{name}\b(?!\s*\()", f"{name}(1)", text, flags=re.IGNORECASE
            )
        text = text.replace("dimension v11(1)(1),v12(1)(1),v22(1)(1)",
                            "dimension v11(1),v12(1),v22(1)")
    if "ldg_r, v11, De, ldg_e, 1)" in text and not re.search(
        r"(?i)dimension\s+v11\s*\(", text
    ):
        text = text.replace(
            "dimension g_r(ldg_r,nt,3),g_e(ldg_e,nt)",
            "dimension g_r(ldg_r,nt,3),g_e(ldg_e,nt)\n"
            "        dimension v11(1),v12(1),v22(1)",
            1,
        )
        for name in ("v11", "v12", "v22"):
            text = re.sub(
                rf"\b{name}\b(?!\s*\()", f"{name}(1)", text, flags=re.IGNORECASE
            )
        text = text.replace("dimension v11(1)(1),v12(1)(1),v22(1)(1)",
                            "dimension v11(1),v12(1),v22(1)")
    if "call pot11(r(i,1),v11,1)" in text and not re.search(
        r"(?i)dimension\s+v11\s*\(", text
    ):
        text = re.sub(
            r"(?im)^(\s*dimension\s+r\(nt,3\),e\(nt\)[^\n]*)$",
            r"\1\n        dimension v11(1),v12(1),v22(1)",
            text,
            count=1,
        )
        for name in ("v11", "v12", "v22"):
            text = re.sub(
                rf"\b{name}\b(?!\s*\()", f"{name}(1)", text, flags=re.IGNORECASE
            )
        text = text.replace("dimension v11(1)(1),v12(1)(1),v22(1)(1)",
                            "dimension v11(1),v12(1),v22(1)")
    if re.search(r"(?im)^\s*ENTRY\s+POT11\s*\(\s*rr\s*,\s*e\s*,\s*n\s*\)", text):
        text = re.sub(r"(?i)\brr\s*\(\s*n\s*,\s*3\s*\)", "rr(1,3)", text)
    text = text.replace(
        "        real g_r(ldg_r, 1, 3)",
        "        double precision g_r(ldg_r, 1, 3)",
    )
    text = text.replace(
        "COMMON /PT2CM/ NDER, NDUM, IDEBUG, NDUM(19)",
        "COMMON /PT2CM/ NDER, NDUM, IDEBUG, IDUM(19)",
    )
    text = re.sub(r"(?im)^[ \t]*integer\s+g_ehfid\s*$", "", text)
    text = re.sub(r"(?im)^[ \t]*data\s+g_ehfid\s*/\s*0\s*/\s*$", "", text)
    if ("CSP(1,1),SCRAP,ANGLE" in text
            and not re.search(r"(?im)^\s*dimension\s+scrap\(", text)):
        text = text.replace(
            "      CHARACTER*20 TITLE",
            "      CHARACTER*20 TITLE\n      DIMENSION SCRAP(22)",
        )
    replacements = {
        "IMPLICIT REAL*8": "IMPLICIT DOUBLE PRECISION",
        "implicit real*8": "implicit double precision",
        "CHARACTER*75": "CHARACTER(75)",
        "CHARACTER*20": "CHARACTER(20)",
        "CHARACTER*3": "CHARACTER(3)",
        "CHARACTER*2": "CHARACTER(2)",
        "CHARACTER*1": "CHARACTER(1)",
        "JC05132003": "",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"(?m)^\*{5,}", lambda match: "!" + match.group(0), text)
    text = re.sub(r"(?m)^[cC](?=\s|$)", "!", text)
    text = re.sub(
        r"(?i)\.\s*(LT|LE|GT|GE|EQ|NE|OR|AND)\s*\.",
        lambda match: f".{match.group(1)}.",
        text,
    )
    text = re.sub(
        r"(?i)(\d)[d]\s+([+-]?)\s*(\d+)",
        lambda match: f"{match.group(1)}D{match.group(2)}{match.group(3)}",
        text,
    )
    text = re.sub(
        r"(\d+\.\d*)\s+(\d{3,})\s+(\d{3,}[dD][+-]?\d+)",
        r"\1\2\3",
        text,
    )
    text = re.sub(
        r"(?im)^([ \t]*do\s*)(\d+)([A-Za-z_][A-Za-z0-9_]*)\s*(?==)",
        r"\1\2 \3",
        text,
    )
    # Some fixed-form files retain historical sequence-field revision marks.
    # They are outside columns 1-72 and are not Fortran source.
    text = re.sub(r"(?m)(^.{0,72})\s+[0-9]+[a-z][0-9]+[a-z][0-9]+\s*$", r"\1", text)
    text = re.sub(
        r"(?im)^(.{6,}?)\s+(?:\d{1,2}[A-Z]{3}\d{2,4}|\d+[A-Z]\d+[A-Z]\d+)\s*$",
        r"\1",
        text,
    )
    lines = []
    for line in text.splitlines(keepends=True):
        body = line.rstrip("\r\n")
        ending = line[len(body):]
        if (
            len(body) > 72
            and body[71].isspace()
            and re.fullmatch(r"[0-9A-Za-z./-]+", body[72:].strip())
        ):
            body = body[:72].rstrip()
        lines.append(body + ending)
    text = "".join(lines)
    return text


def strip_inline_comment(text: str) -> str:
    quote = ""
    index = 0
    while index < len(text):
        character = text[index]
        if quote:
            if character == quote:
                if index + 1 < len(text) and text[index + 1] == quote:
                    index += 2
                    continue
                quote = ""
        elif character in {"'", '"'}:
            quote = character
        elif character == "!":
            return text[:index].rstrip()
        index += 1
    return text.rstrip()


def fixed_to_free(text: str) -> str:
    """Convert normalized fixed-form source to free form before transformation."""
    output: list[str] = []
    statement: str | None = None
    pending_comments: list[str] = []

    def flush() -> None:
        nonlocal statement
        if statement is not None:
            output.append(statement.rstrip())
            statement = None
        if pending_comments:
            output.extend(pending_comments)
            pending_comments.clear()

    for raw in text.splitlines():
        raw = raw.expandtabs(8)
        if not raw.strip():
            if statement is None:
                output.append("")
            else:
                pending_comments.append("")
            continue
        if raw[0] in "cCdD*!":
            if statement is None:
                output.append("!" + raw[1:] if raw[0] != "!" else raw)
            else:
                pending_comments.append("!" + raw[1:] if raw[0] != "!" else raw)
            continue
        raw = strip_inline_comment(raw)
        if not raw:
            continue
        if re.match(
            r"^ {5}(?:IMPLICIT|INTEGER|REAL|DOUBLE|CHARACTER|LOGICAL|"
            r"DIMENSION|PARAMETER|COMMON|SAVE|DATA|CALL|RETURN|ENTRY|"
            r"SUBROUTINE|FUNCTION|BLOCK|DO|IF|GOTO|GO\s+TO)\b",
            raw,
            flags=re.I,
        ):
            raw = " " + raw
        padded = raw + " " * max(0, 6 - len(raw))
        continuation = padded[5] not in {" ", "0"}
        body = padded[6:] if len(padded) > 6 else ""
        if continuation:
            if statement is None:
                raise RuntimeError("fixed-form continuation has no leading statement")
            statement += body.rstrip()
            continue
        flush()
        label = padded[:5].strip()
        statement = (label + " " if label else "") + body.rstrip()
    flush()
    return "\n".join(output) + "\n"


def join_free_continuations(text: str) -> str:
    output: list[str] = []
    statement = ""
    continuing = False
    for raw in text.splitlines():
        stripped = raw.strip()
        if continuing and (not stripped or stripped.startswith("!")):
            continue
        if not continuing and (not stripped or stripped.startswith("!")):
            output.append(raw)
            continue
        content = raw.lstrip()
        if continuing and content.startswith("&"):
            content = content[1:].lstrip()
        has_continuation = content.rstrip().endswith("&")
        if has_continuation:
            content = content.rstrip()[:-1]
        statement += content
        continuing = has_continuation
        if not continuing:
            output.append(statement)
            statement = ""
    if statement:
        output.append(statement)
    return "\n".join(output) + "\n"


QUIET = r"""
C potlibc generated quiet, fixed-unit POTLIB initialization.
      SUBROUTINE POTLIBC_POTINFO_QUIET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NATOM=25,ISURF=5,JSURF=15)
      CHARACTER(75) REF(5)
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     + IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,NULBL(NATOM),
     + NFLAG(20),NASURF(ISURF+1,ISURF+1),NDER
      COMMON/UTILCM/ DGSCARTNU(NATOM,3),DESCARTNU(NATOM,3,ISURF),
     + DIJCARTNU(NATOM,3,JSURF),CNVRTD,CNVRTE,CNVRTDE,
     + IREORDER,KSDIAG,KEDIAG,KSOFFD,KEOFFD
      NDER=MDER
      KSDIAG=0
      KEDIAG=0
      DO I=2,ISURF+1
        IF(NASURF(I,I).NE.0) THEN
          KEDIAG=I-1
          IF(KSDIAG.EQ.0) KSDIAG=I-1
        ENDIF
      ENDDO
      KSOFFD=0
      KEOFFD=0
      K=0
      DO I=1,ISURF
        DO J=I+1,ISURF+1
          K=K+1
          IF(NASURF(I,J)+NASURF(J,I).NE.0) THEN
            KEOFFD=K
            IF(KSOFFD.EQ.0) KSOFFD=K
          ENDIF
        ENDDO
      ENDDO
      CNVRTD=1.0D0
      CNVRTE=1.0D0
      CNVRTDE=1.0D0
      IREORDER=0
      RETURN
      END
"""


ADAPTER = r"""
subroutine potlib_plugin_initialize() bind(c)
  implicit none
  call prepot()
__DUAL_INITIALIZE__
end subroutine

subroutine potlib_plugin_metadata(natoms, derivatives, atomic_numbers, reference) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), intent(out) :: natoms, derivatives, atomic_numbers(25)
  real(c_double), intent(out) :: reference(75)
  integer, parameter :: max_atoms=25
  character(75) :: ref(5)
  double precision :: cartnu(max_atoms,3)
  integer :: indexes(max_atoms), irctnt, native_atoms, icartr, mder, msurf, i, j, k
  common /infocm/ cartnu,indexes,irctnt,native_atoms,icartr,mder,msurf,ref
  natoms=native_atoms
  derivatives=mder
  do i=1,max_atoms
    atomic_numbers(i)=indexes(i)
    do j=1,3
      k=(i-1)*3+j
      reference(k)=cartnu(i,j)
    end do
  end do
end subroutine

subroutine potlib_plugin_evaluate(n, x, energy, force) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none
__DUAL_INTERFACE__
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*)
  real(c_double), intent(out) :: energy(*), force(*)
  integer, parameter :: max_atoms=25, isurf=5, jsurf=15
  double precision :: pengygs,pengyes(isurf),pengyij(jsurf)
  double precision :: dgscart(max_atoms,3),descart(max_atoms,3,isurf)
  double precision :: dijcart(max_atoms,3,jsurf),cart(max_atoms,3),anuzero
  integer :: nullbl(max_atoms),nflag(20),nasurf(isurf+1,isurf+1),nder
  character(75) :: ref(5)
  double precision :: cartnu(max_atoms,3)
  integer :: indexes(max_atoms),irctnt,natoms,icartr,mder,msurf
  integer :: sample,atom,axis,k
  common /usrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /usricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  common /infocm/ cartnu,indexes,irctnt,natoms,icartr,mder,msurf,ref
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+(atom-1)*3+axis
        cart(atom,axis)=x(k)
      end do
    end do
    call pot()
    energy(sample)=pengygs
    if (mder.ge.1) then
      do atom=1,natoms
        do axis=1,3
          k=(sample-1)*3*natoms+(atom-1)*3+axis
          force(k)=-dgscart(atom,axis)
        end do
      end do
    endif
  end do
__FORCE_DUAL__
end subroutine
"""

MODULE_ADAPTER = r"""
subroutine potlib_plugin_initialize() bind(c)
  implicit none
end subroutine

subroutine potlib_plugin_metadata(natoms, derivatives, atomic_numbers, reference) bind(c)
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), intent(out) :: natoms, derivatives, atomic_numbers(25)
  real(c_double), intent(out) :: reference(75)
  integer :: i
  natoms=__NATOMS__
  derivatives=1
  atomic_numbers=0
  atomic_numbers(1:__NATOMS__)=[__ATOMIC_NUMBERS__]
  do i=1,75
    reference(i)=0.0_c_double
  end do
end subroutine

subroutine potlib_plugin_evaluate(n, x, energy, force) bind(c)
  use iso_c_binding, only: c_double, c_int
  use __MODULE__, only: pot
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*)
  real(c_double), intent(out) :: energy(*), force(*)
  real(c_double) :: xc(__NATOMS__), yc(__NATOMS__), zc(__NATOMS__)
  real(c_double) :: dx(__NATOMS__), dy(__NATOMS__), dz(__NATOMS__)
  integer :: sample, atom, base
  do sample=1,n
    base=(sample-1)*__NDIM__
    do atom=1,__NATOMS__
      xc(atom)=x(base+(atom-1)*3+1)
      yc(atom)=x(base+(atom-1)*3+2)
      zc(atom)=x(base+(atom-1)*3+3)
    end do
    call pot(xc,yc,zc,energy(sample),dx,dy,dz)
    do atom=1,__NATOMS__
      force(base+(atom-1)*3+1)=-dx(atom)
      force(base+(atom-1)*3+2)=-dy(atom)
      force(base+(atom-1)*3+3)=-dz(atom)
    end do
  end do
end subroutine
"""

DUAL_INTERFACE = r"""
  interface
    subroutine potlib_plugin_dual_force(n, x, force) bind(c)
      import c_double, c_int
      integer(c_int), value :: n
      real(c_double), intent(in) :: x(*)
      real(c_double), intent(out) :: force(*)
    end subroutine
  end interface
"""

DUAL_ADAPTER = r"""
subroutine potlib_plugin_dual_force(n, x, force) bind(c)
  use dual_mod, only: dual
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*)
  real(c_double), intent(out) :: force(*)
  integer, parameter :: max_atoms=25,isurf=5,jsurf=15,natoms=__NATOMS__
  type(dual) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual) :: dgscart(max_atoms,3),descart(max_atoms,3,isurf)
  type(dual) :: dijcart(max_atoms,3,jsurf),cart(max_atoms,3),anuzero
  integer :: nullbl(max_atoms),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: sample,atom,axis,k,target
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do sample=1,n
    do target=1,3*natoms
      do atom=1,natoms
        do axis=1,3
          k=(sample-1)*3*natoms+(atom-1)*3+axis
          cart(atom,axis)=dual(x(k),0.0d0)
        end do
      end do
      atom=(target-1)/3+1
      axis=mod(target-1,3)+1
      cart(atom,axis)%du=1.0d0
      call zpot()
      k=(sample-1)*3*natoms+target
      force(k)=-pengygs%du
    end do
  end do
end subroutine

subroutine potlib_plugin_energy_directional(n, x, direction, derivative) bind(c)
  use dual_mod, only: dual
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*), direction(*)
  real(c_double), intent(out) :: derivative(*)
  integer, parameter :: max_atoms=25,isurf=5,jsurf=15,natoms=__NATOMS__
  type(dual) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual) :: dgscart(max_atoms,3),descart(max_atoms,3,isurf)
  type(dual) :: dijcart(max_atoms,3,jsurf),cart(max_atoms,3),anuzero
  integer :: nullbl(max_atoms),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: sample,atom,axis,k
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+(atom-1)*3+axis
        cart(atom,axis)=dual(x(k),direction(k))
      end do
    end do
    call zpot()
    derivative(sample)=pengygs%du
  end do
end subroutine

subroutine potlib_plugin_hvp(n, x, direction, product) bind(c)
  use dual_mod, only: dual
  use iso_c_binding, only: c_double, c_int
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*), direction(*)
  real(c_double), intent(out) :: product(*)
  integer, parameter :: max_atoms=25,isurf=5,jsurf=15,natoms=__NATOMS__
  type(dual) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual) :: dgscart(max_atoms,3),descart(max_atoms,3,isurf)
  type(dual) :: dijcart(max_atoms,3,jsurf),cart(max_atoms,3),anuzero
  integer :: nullbl(max_atoms),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: sample,atom,axis,k
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+(atom-1)*3+axis
        cart(atom,axis)=dual(x(k),direction(k))
      end do
    end do
    call zpot()
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+(atom-1)*3+axis
        product(k)=dgscart(atom,axis)%du
      end do
    end do
  end do
end subroutine
"""


def adapted_source(text: str, fixed_unit_init: bool, cbe: bool) -> str:
    text = normalize(text)
    if fixed_unit_init:
        text, potinfo_count = re.subn(
            ACTIVE_CALL.format("potinfo"),
            lambda m: re.sub(
                "(?i)potinfo", "POTLIBC_POTINFO_QUIET", m.group(1)
            ),
            text,
        )
        text, ancvrt_count = re.subn(
            ACTIVE_CALL.format("ancvrt"),
            lambda m: "C potlibc replaced ANCVRT with fixed-unit initialization",
            text,
        )
        if potinfo_count != 1 or ancvrt_count != 1:
            raise RuntimeError("fixed-unit initialization metadata changed")
    if cbe:
        text, count = re.subn(
            r"(?i)DATA\s+ANUZERO\s*/\s*0\.0D0\s*/",
            "DATA ANUZERO /-0.7098324298858642578125D0/",
            text,
        )
        if count != 1:
            raise RuntimeError("expected one CBE ANUZERO DATA initializer")
    return text + (QUIET if fixed_unit_init else "")


def harvest(text: str) -> tuple[set[str], set[str]]:
    active = "".join(
        line
        for line in text.splitlines(keepends=True)
        if not line.lstrip().startswith(("!", "*"))
    )
    identifiers = set(
        re.findall(r"(?i)^\s*subroutine\s+([A-Za-z_][\w]*)", active, re.M)
    )
    identifiers |= set(
        re.findall(
            r"(?i)^\s*(?:(?:double\s+precision|real\*8|complex\*16)\s+)?"
            r"function\s+([A-Za-z_][\w]*)",
            active,
            re.M,
        )
    )
    identifiers |= set(re.findall(r"(?i)\bcall\s+([A-Za-z_][\w]*)", active))
    commons = set(
        re.findall(r"(?i)common\s*/\s*([A-Za-z_][\w]*)\s*/", active)
    )
    return identifiers, commons


def dual_literals(text: str) -> str:
    lines = text.splitlines(keepends=True)
    return "".join(
        line if "parameter" in line.lower() else dual_literals_line(line)
        for line in lines
    )


def dual_literals_line(text: str) -> str:
    pattern = re.compile(r"(?<![A-Za-z])((?:\d+\.\d*|\.\d+)[dD][+-]?\d+)")
    parts: list[str] = []
    last = 0
    for match in pattern.finditer(text):
        parts.append(text[last : match.start()])
        parts.append(f"({match.group(1)},0d0)")
        last = match.end()
    parts.append(text[last:])
    text = "".join(parts)
    text = re.sub(r"([,/])\s*-\(([^,\n]+),0d0\)", r"\1(-\2,0d0)", text)
    text = re.sub(r"([,/])\s*\+\(([^,\n]+),0d0\)", r"\1(+\2,0d0)", text)
    return text


def to_dual(text: str, identifiers: set[str], commons: set[str]) -> str:
    if re.search(r"(?im)^[ \t]*sign\s*=", text):
        text = re.sub(r"(?i)\bsign\b", "potlibc_sign_value", text)
        text = re.sub(
            r"(?i)\bpotlibc_sign_value\s*\(",
            "sign(",
            text,
        )
    text = re.sub(
        r"(?im)^([ \t]*)(?:(\d+)\s+)?IF\s*\((.*)\)\s*"
        r"(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*$",
        lambda match: (
            f"{match.group(1)}{(match.group(2) + ' ') if match.group(2) else ''}"
            f"IF (({match.group(3)}) .LT. 0) GOTO {match.group(4)}\n"
            f"{match.group(1)}IF (({match.group(3)}) .EQ. 0) GOTO {match.group(5)}\n"
            f"{match.group(1)}GOTO {match.group(6)}"
        ),
        text,
    )
    text = re.sub(
        r"(?i)\bBLOCK\s+DATA\s+PTPACM\b", "BLOCK DATA zptpacm", text
    )
    for name in sorted(identifiers, key=len, reverse=True):
        text = re.sub(rf"\b{name}\b", "z" + name.lower(), text, flags=re.I)
    for name in sorted(commons, key=len, reverse=True):
        text = re.sub(rf"/\s*{name}\s*/", f"/z{name.lower()}/", text, flags=re.I)
    text = re.sub(
        r"(?im)^([ \t]*module\s+(?!procedure\b)[A-Za-z_][A-Za-z0-9_]*\s*\n)",
        r"\1      use dual_mod\n",
        text,
    )
    text = re.sub(
        r"(?im)^([ \t]*(?:(?:double\s+precision|real\*8|complex\*16)\s+)?"
        r"(?:subroutine|function|block\s+data)\b[^\n]*"
        r"(?:\n[ \t]{5}[^ 0\n][^\n]*)*\n)",
        r"\1      use dual_mod\n",
        text,
    )
    text = re.sub(
        r"(?i)\bdouble\s+precision\s*,\s*parameter\b",
        "REAL*8, parameter",
        text,
    )
    text = re.sub(
        r"(?im)^([ \t]*)double\s+precision([^\n]*=[^\n]*)$",
        r"\1REAL*8\2",
        text,
    )

    def split_declarators(body: str) -> list[str]:
        declarations: list[str] = []
        start = 0
        depth = 0
        for index, character in enumerate(body):
            if character in "([":
                depth += 1
            elif character in ")]":
                depth -= 1
            elif character == "," and depth == 0:
                declarations.append(body[start:index].strip())
                start = index + 1
        declarations.append(body[start:].strip())
        return declarations

    def rewrite_real_declaration(match: re.Match[str]) -> str:
        indent, body = match.group(1), match.group(2)
        declarations = split_declarators(body)
        unit_starts = [
            candidate.start()
            for candidate in re.finditer(
                r"(?im)^[ \t]*(?:module\s+(?!procedure\b)|subroutine|"
                r"(?:[a-z][a-z0-9_*( ),]+\s+)?function|program|block\s+data)\b",
                text,
            )
        ]
        start = max(
            (position for position in unit_starts if position <= match.start()),
            default=0,
        )
        end = min(
            (position for position in unit_starts if position > match.start()),
            default=len(text),
        )
        local_unit = text[start:end]
        local_parameters = {
            name.lower()
            for line in local_unit.splitlines()
            if re.search(r"(?i)\bparameter\s*(?:\(|::)", line)
            and not re.match(r"^\s*!", line)
            for name in re.findall(
                r"([A-Za-z_][A-Za-z0-9_]*)\s*=", line, flags=re.I
            )
        }
        constants: list[str] = []
        integers: list[str] = []
        variables: list[str] = []
        for declaration in declarations:
            name_match = re.match(
                r"(?:\s*::\s*)?([A-Za-z_][A-Za-z0-9_]*)", declaration
            )
            if (
                name_match
                and (
                    name_match.group(1).lower() in local_parameters
                    or "=" in declaration
                )
            ):
                constants.append(declaration)
            elif name_match and re.fullmatch(
                r"ind\d+", name_match.group(1), flags=re.I
            ):
                integers.append(declaration)
            else:
                variables.append(declaration)
        lines = []
        if constants:
            lines.append(
                indent
                + "REAL*8 :: "
                + ",".join(
                    re.sub(r"^\s*::\s*", "", declaration)
                    for declaration in constants
                )
            )
        if integers:
            lines.append(indent + "INTEGER " + ",".join(integers))
        if variables:
            lines.append(indent + "type(dual) " + ",".join(variables))
        return "\n".join(lines)

    text = re.sub(
        r"(?im)^([ \t]*)REAL\*8\s+(?!FUNCTION\b)([^\n]+)$",
        rewrite_real_declaration,
        text,
    )
    text = re.sub(
        r"(?im)^([ \t]*)REAL\*8(\s*,(?![^\n]*\bparameter\b)(?![^\n]*=)[^\n]+)$",
        r"\1type(dual)\2",
        text,
    )
    text = re.sub(
        r"(?im)^([ \t]*)REAL\*8\s+FUNCTION\b",
        r"\1type(dual) FUNCTION",
        text,
    )
    text = re.sub(
        r"(?im)^([ \t]*)DOUBLE\s+PRECISION\s+(?!FUNCTION\b)([^\n]+)$",
        rewrite_real_declaration,
        text,
    )
    text = re.sub(
        r"implicit\s+double\s+precision", "implicit complex*16", text, flags=re.I
    )
    text = re.sub(r"\bdouble\s+precision\b", "complex*16", text, flags=re.I)
    text = re.sub(
        r"(?im)^([ \t]*)implicit\s+complex\*16\s*(\([^)\n]*\))",
        r"\1implicit type(dual) \2",
        text,
    )
    text = text.replace("complex*16", "type(dual)")
    text = re.sub(r"\(([^,\n()]+),0d0\)", r"dual(\1,0d0)", text)
    text = re.sub(
        r"([,/])\s*-\s*dual\(([^,\n]+),0d0\)",
        r"\1dual(-\2,0d0)",
        text,
    )
    data_statement = False
    data_slashes = 0
    data_lines: list[str] = []
    for line in text.splitlines(keepends=True):
        starts_data = bool(re.match(r"(?i)^(?![c*!])[^\n]*\bDATA\b", line))
        if starts_data:
            data_statement = True
            data_slashes = 0
        if data_statement:
            line = re.sub(
                r"-\s*dual\(([^,\n]+),0d0\)",
                r"dual(-\1,0d0)",
                line,
            )
            data_slashes += line.count("/")
            if data_slashes >= 2:
                data_statement = False
        data_lines.append(line)
    text = "".join(data_lines)
    text = re.sub(
        r"(?i)\bD(MAX1|MIN1)\b",
        lambda match: match.group(1)[:-1],
        text,
    )
    text = re.sub(
        r"(?i)\bD(SQRT|EXP|LOG10|LOG|ABS|SINH?|COSH?|TANH?|ASIN|ACOS|ATAN|ERFC?|MOD)\b",
        r"\1",
        text,
    )
    text = re.sub(r"(?i)\bIDINT\b", "INT", text)
    text = re.sub(r"(?i)\bIDNINT\b", "NINT", text)
    text = re.sub(r"(?im)^[ \t]*INTRINSIC\s+DBLE\s*$", "", text)

    explicitly_declared = {
        name.lower()
        for declaration in re.findall(
            r"(?im)^[ \t]*(?:INTEGER|LOGICAL|CHARACTER(?:\([^)]*\)|\*\d+)?|"
            r"REAL(?:\*\d+)?|DOUBLE\s+PRECISION|COMPLEX(?:\*\d+)?|"
            r"type\s*\(\s*dual\s*\))\s*(?:::\s*)?([^\n]+)$",
            text,
        )
        for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", declaration)
    }
    def declare_real_parameters(match: re.Match[str]) -> str:
        line = match.group(0)
        unit_starts = [
            candidate.start()
            for candidate in re.finditer(
                r"(?im)^[ \t]*(?:subroutine|(?:[a-z][a-z0-9_*( ),]+\s+)?"
                r"function|program|block\s+data)\b",
                text[:match.start()],
            )
        ]
        local_text = text[unit_starts[-1] if unit_starts else 0:match.start()]
        locally_declared = {
            name.lower()
            for declaration in re.findall(
                r"(?im)^[ \t]*(?:INTEGER|LOGICAL|CHARACTER(?:\([^)]*\)|\*\d+)?|"
                r"REAL(?:\*\d+)?|DOUBLE\s+PRECISION|COMPLEX(?:\*\d+)?)"
                r"\s*(?:::\s*)?([^\n]+)$",
                local_text,
            )
            for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", declaration)
        }
        names = [
            name
            for name in re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
            if (name[0].lower() < "i" or name[0].lower() > "n")
            and name.lower() not in locally_declared
        ]
        if not names:
            return line
        return "      REAL*8 " + ",".join(names) + "\n" + line

    text = re.sub(
        r"(?im)^(?![c*!])[^\n]*\bPARAMETER\s*\([^\n]*$",
        declare_real_parameters,
        text,
    )
    text = re.sub(
        r"PARAMETER\(\s*CMU_0\s*=\s*dual\([^)]*\)\*PI\*dual\([^)]*\)\)",
        "PARAMETER(CMU_0=dual(1.2566370614359173d-06,0d0))",
        text,
    )
    text = re.sub(
        r"PARAMETER\(CEPSILON_0\s*=\s*dual\([^)]*\)/\(CMU_0\*CLIGHT\*\*2\)\)",
        "PARAMETER(CEPSILON_0=dual(8.8541878176203899d-12,0d0))",
        text,
    )
    def data_line(match: re.Match[str]) -> str:
        body = match.group(0)
        targets = body[re.search(r"(?i)\bDATA\b", body).end():].split("/", 1)[0]
        target_names = {
            name.lower()
            for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", targets)
        }
        unit_starts = [
            candidate.start()
            for candidate in re.finditer(
                r"(?im)^[ \t]*(?:subroutine|(?:[a-z][a-z0-9_*( ),]+\s+)?"
                r"function|program|block\s+data)\b",
                text[:match.start()],
            )
        ]
        local_text = text[unit_starts[-1] if unit_starts else 0:match.start()]
        local_real = {
            name.lower()
            for declaration in re.findall(
                r"(?im)^[ \t]*REAL(?:\*\d+)?\s*(?:,[^:]*)?"
                r"(?:::\s*)?([^\n]+)$",
                local_text,
            )
            for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", declaration)
        }
        if target_names & local_real:
            return body
        def wrap(number: re.Match[str]) -> str:
            token = number.group(0)
            if "dual(" in body[max(0, number.start() - 5) : number.start()]:
                return token
            suffix = "" if "d" in token.lower() else "d0"
            return f"dual({token}{suffix},0d0)"
        return re.sub(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[dD][+-]?\d+)?", wrap, body)
    text = re.sub(
        r"(?im)^[^\n]*\bDATA\b[^\n]*",
        data_line,
        text,
    )
    return text


def to_dual2(text: str, identifiers: set[str], commons: set[str]) -> str:
    transformed = to_dual(text, identifiers, commons)
    transformed = transformed.replace("use dual_mod", "use dual2_mod")
    transformed = transformed.replace("type(dual)", "type(dual2)")
    transformed = re.sub(r"\bdual\(", "dual2(", transformed)
    return re.sub(
        r"\bdual2\(([^,()\n]+),\s*0d0\)",
        r"dual2(\1,0d0,0d0,0d0)",
        transformed,
    )


_JET_LITERAL = re.compile(
    r"(?<!jet3\()(?<!jet3\(-)(?<!jet3\(\+)(?<![\w.'])"
    r"([+-]?)((?:\d+\.\d*|\.\d+)[dD][+-]?\d+)(?![\w.'])"
)


def _wrap_real_literals(text: str) -> str:
    """Make every double-precision literal a jet constant.

    The operation count of a jet is the count of jet operations, and a
    constant operand is one of them, as in the reference CBE kernel that this
    transform reproduces.
    """
    out = []
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped[:1] in ("!", "c", "C", "*") and not stripped[:2].lower() in (
            "co", "ca", "ch",
        ):
            out.append(line)
            continue
        if re.match(r"(?i)parameter\b", stripped):
            out.append(line)
            continue
        code, bang, comment = line.partition("!")
        code = _JET_LITERAL.sub(
            lambda m: f"{m.group(1)}jet3({m.group(2)},0d0,0d0,0d0)", code
        )
        out.append(code + bang + comment)
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def to_jet3(text: str, identifiers: set[str], commons: set[str]) -> str:
    """Retype a to_dual2 source on third-order jets under j-prefixed names."""
    for name in sorted(identifiers, key=len, reverse=True):
        text = re.sub(
            rf"\bz{re.escape(name.lower())}\b",
            "j" + name.lower(),
            text,
            flags=re.I,
        )
    for name in sorted(commons, key=len, reverse=True):
        text = re.sub(
            rf"/\s*z{re.escape(name.lower())}\s*/",
            f"/j{name.lower()}/",
            text,
            flags=re.I,
        )
    text = re.sub(r"(?i)\bBLOCK\s+DATA\s+zptpacm\b", "BLOCK DATA jptpacm", text)
    text = text.replace("use dual2_mod", "use jet3_mod")
    text = text.replace("type(dual2)", "type(jet3)")
    text = re.sub(r"\bdual2\(", "jet3(", text)
    return _wrap_real_literals(text)


def emit(args: argparse.Namespace) -> None:
    vendor = Path(args.vendor)
    entries = inventory(vendor)
    entry = next(
        item for item in entries if item.get("wrapped") and item.get("id") == args.surface
    )
    source = vendor / str(entry["path"])
    cbe = str(entry["path"]) == "src/CH5/ch52008.f"
    fixed_unit_init = bool(entry.get("fixed_unit_init"))
    Path(args.output_source).write_text(
        adapted_source(source.read_text(encoding="latin-1"), fixed_unit_init, cbe)
    )
    Path(args.output_utility).write_text(normalize((vendor / "utilities/utility.f").read_text()))
    dual = bool(entry.get("hvp"))
    module_match = re.search(
        r"(?im)^\s*module\s+([a-z][a-z0-9_]*)\s*$",
        source.read_text(encoding="latin-1"),
    )
    if module_match:
        module_atoms = 4 if "/N2O2/" in f"/{entry['path']}" else 3
        atomic_numbers = "7,7,8,8" if module_atoms == 4 else "7,7,8"
        adapter = MODULE_ADAPTER.replace("__MODULE__", module_match.group(1))
        adapter = adapter.replace("__NATOMS__", str(module_atoms))
        adapter = adapter.replace("__NDIM__", str(3 * module_atoms))
        adapter = adapter.replace("__ATOMIC_NUMBERS__", atomic_numbers)
    else:
        adapter = ADAPTER.replace(
            "__DUAL_INITIALIZE__", "  call zprepot()" if dual else ""
        )
        adapter = adapter.replace(
            "__DUAL_INTERFACE__", DUAL_INTERFACE if dual else ""
        ).replace(
            "__FORCE_DUAL__",
            "  if (mder.lt.1) call potlib_plugin_dual_force(n,x,force)"
            if dual
            else "",
        )
    if dual:
        adapter += DUAL_ADAPTER.replace("__NATOMS__", str(entry["atoms"]))
    Path(args.output_adapter).write_text(adapter)
    if args.dual_source is not None:
        if not dual:
            raise RuntimeError("dual outputs requested for a surface without HVP support")
        real_source = adapted_source(source.read_text(encoding="latin-1"), fixed_unit_init, cbe)
        real_utility = normalize((vendor / "utilities/utility.f").read_text())
        source_ids, source_commons = harvest(real_source)
        utility_ids, utility_commons = harvest(real_utility)
        identifiers = source_ids | utility_ids
        commons = source_commons | utility_commons
        Path(args.dual_source).write_text(to_dual(real_source, identifiers, commons))
        Path(args.dual_utility).write_text(to_dual(real_utility, identifiers, commons))


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    inv = sub.add_parser("inventory")
    inv.add_argument("vendor")
    inv.add_argument("output")
    inv.add_argument("--meson", action="store_true")
    inv.add_argument("--notice", action="store_true")
    inv.add_argument("--provenance", action="store_true")
    generate = sub.add_parser("emit")
    generate.add_argument("vendor")
    generate.add_argument("surface")
    generate.add_argument("output_source")
    generate.add_argument("output_utility")
    generate.add_argument("output_adapter")
    generate.add_argument("--dual-source")
    generate.add_argument("--dual-utility")
    args = parser.parse_args()
    if args.command == "inventory":
        data = inventory(Path(args.vendor))
        if args.meson:
            print(
                "\n".join(
                    f"{item['id']}|{item['path']}|{int(bool(item['hvp']))}"
                    for item in data
                    if item["wrapped"]
                )
            )
        elif args.notice:
            lines = [
                "potlibc",
                "Copyright 2026 Rohit Goswami",
                "",
                "This distribution includes potential-energy surfaces from POTLIB:",
                "",
                "Y. Shu, Z. Varga, A. W. Jasper, B. C. Garrett, J. Espinosa-García,",
                "J. C. Corchado, R. J. Duchovic, M. Volobuev, B. J. Lynch, K. Yang,",
                "T. C. Allison, A. F. Wagner, and D. G. Truhlar, POTLIB,",
                "https://comp.chem.umn.edu/potlib",
                "",
                "Y. Shu, Z. Varga, A. W. Jasper, J. Espinosa-García, J. C. Corchado,",
                "and D. G. Truhlar, “PotLib 2023: New version of a potential energy",
                "surface library for chemical systems,” Mendeley Data, version 1 (2023).",
                "https://doi.org/10.17632/chxxhnnt4p.1",
                "",
                "Included surface file-header references:",
            ]
            for item in data:
                source_fields = fields(
                    (Path(args.vendor) / str(item["path"])).read_text(
                        errors="replace"
                    )
                )
                reference = source_fields.get(
                    "reference",
                    source_fields.get(
                        "references",
                        "(file header has no Reference(s) field; "
                        "see the cited vendor path)",
                    ),
                )
                lines.append(f"- {item['path']}: {reference}")
            lines += [
                "",
                "utilities/utility.f was written by R. J. Duchovic and A. F. Wagner",
                "and includes corrections by J. Corchado.",
            ]
            Path(args.output).write_text("\n".join(lines) + "\n")
        elif args.provenance:
            vendor_path = Path(args.vendor)
            files = sorted(path for path in vendor_path.rglob("*") if path.is_file())
            lines = [
                "<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->",
                "<!-- SPDX-License-Identifier: BSD-3-Clause -->",
                "# Provenance",
                "",
                "The vendor tree is a byte-for-byte copy of the PotLib 2023 archive:",
                "",
                "> Y. Shu, Z. Varga, A. Jasper, J. Espinosa-García, J. C. Corchado,",
                "> D. G. Truhlar, “PotLib 2023: New version of a potential energy",
                "> surface library for chemical systems,” Mendeley Data, version 1",
                "> (2023), <https://doi.org/10.17632/chxxhnnt4p.1>.",
                "",
                "| Repository path | SHA-256 |",
                "| --- | --- |",
            ]
            for path in files:
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                lines.append(
                    f"| `upstream/potlib2023/{path.relative_to(vendor_path).as_posix()}` "
                    f"| `{digest}` |"
                )
            Path(args.output).write_text("\n".join(lines) + "\n")
        else:
            output = Path(args.output)
            if output.suffix == ".md":
                wrapped = sum(bool(item["wrapped"]) for item in data)
                lines = [
                    "<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->",
                    "<!-- SPDX-License-Identifier: BSD-3-Clause -->",
                    "# Surface inventory",
                    "",
                    f"Generated from POTLIB headers: {wrapped} wrapped, "
                    f"{len(data) - wrapped} unwrapped, {len(data)} total Fortran entries.",
                    "",
                    "| Vendor path | Surface id | Interface | Status | Reason |",
                    "| --- | --- | --- | --- | --- |",
                ]
                for item in data:
                    lines.append(
                        f"| `{item['path']}` | `{item.get('id', '')}` | "
                        f"{item['interface']} | "
                        f"{'wrapped' if item['wrapped'] else 'unwrapped'} | "
                        f"{item['reason']} |"
                    )
                output.write_text("\n".join(lines) + "\n")
            else:
                output.write_text(json.dumps(data, indent=2) + "\n")
    else:
        emit(args)


if __name__ == "__main__":
    main()
