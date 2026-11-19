#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Generate a dual2 POTLIB source, utility, and plugin adapter."""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True

try:
    from .catalog import (
        adapted_source,
        fixed_to_free,
        harvest,
        join_free_continuations,
        normalize,
        to_dual2,
    )
    from .surface_inventory import inventory
except ImportError:
    from catalog import (  # type: ignore[no-redef]
        adapted_source,
        fixed_to_free,
        harvest,
        join_free_continuations,
        normalize,
        to_dual2,
    )
    from surface_inventory import inventory  # type: ignore[no-redef]


PAIR_FAMILIES = {"3-XS", "3-1S"}
DISTANCE_FAMILIES = {"3-1V", "3-2V"}
POTLIB_FAMILY = "POTLIB-2001"
FIXED_SUFFIXES = {".f", ".F"}
ACTIVE_CALL = r"(?im)^(?![c*!])([ \t]*call\s+{}\s*(?:\([^)]*\))?[^\n]*)$"


@dataclass(frozen=True)
class SurfacePlan:
    record: dict[str, Any]
    source: Path
    atoms: int | None
    variable_atoms: bool
    strategy: str
    entry_name: str
    entry_arguments: tuple[str, ...]
    module_name: str | None
    common_name: str | None = None


def _active_call_count(text: str, name: str) -> int:
    return len(re.findall(ACTIVE_CALL.format(re.escape(name)), text))


def _source_text(path: Path) -> str:
    return path.read_text(encoding="latin-1")


def _free_source(text: str, suffix: str, *, normalize_source: bool = True) -> str:
    normalized = normalize(text) if normalize_source else text
    # catalog.normalize() splits this identifier in brh2secd.f.
    normalized = re.sub(r"(?i)\bdo4\s+z(?=\s*\()", "do4z", normalized)
    normalized = re.sub(
        r"(?i)(?<![a-z0-9_])(\d+\.\d*|\.\d+)\s*e\s*([+-]?\d+)",
        r"\1D\2",
        normalized,
    )
    if suffix in FIXED_SUFFIXES:
        normalized = fixed_to_free(normalized)
    return join_free_continuations(normalized)


def _module_name(text: str) -> str | None:
    match = re.search(
        r"(?im)^\s*module\s+(?!procedure\b)([a-z][a-z0-9_]*)\s*$", text
    )
    return match.group(1) if match else None


def _entry_module(text: str, entry_name: str) -> str | None:
    current: str | None = None
    for line in text.splitlines():
        start = re.match(
            r"(?i)^\s*module\s+(?!procedure\b)([a-z][a-z0-9_]*)\s*$", line
        )
        if start:
            current = start.group(1)
            continue
        if re.match(r"(?i)^\s*end\s+module\b", line):
            current = None
            continue
        if current and re.match(
            rf"(?i)^\s*(?:(?:double\s+precision|real\*8)\s+)?"
            rf"(?:subroutine|function)\s+{re.escape(entry_name)}\b",
            line,
        ):
            return current
    return None


def _implementation_text(record: dict[str, Any], vendor: Path) -> str:
    archive = record.get("archive_source")
    if archive is None:
        return _source_text(vendor / record["path"])
    with tarfile.open(vendor / archive["path"]) as package:
        member = package.extractfile(archive["member"])
        if member is None:
            raise ValueError(
                f"{record['id']}: {archive['path']} has no {archive['member']}"
            )
        text = member.read().decode("latin-1")
    first_routine = re.search(r"(?im)^\s*subroutine\s+", text)
    if first_routine is None:
        raise ValueError(f"{record['id']}: archived implementation has no subroutine")
    return text[first_routine.start():]


def _embed_parameter_files(
    text: str, record: dict[str, Any], vendor: Path
) -> str:
    helpers: list[str] = []
    for number, reference in enumerate(record["external_data_files"], 1):
        if not reference.get("literal"):
            continue
        name = reference["name"]
        pattern = re.compile(
            rf"(?im)^(\s*)open\s*\(([^)\n]*\bfile\s*=\s*"
            rf"(['\"]){re.escape(name)}\3[^)\n]*)\)\s*$"
        )
        matches = list(pattern.finditer(text))
        if not matches:
            continue
        contents: list[str] | None = None
        if reference["resolved_paths"]:
            data_path = vendor / reference["resolved_paths"][0]
            if not tarfile.is_tarfile(data_path):
                contents = data_path.read_text(encoding="latin-1").splitlines()
        else:
            source_directory = (vendor / record["path"]).parent
            alternatives = sorted(source_directory.glob(name + ".*"))
            if alternatives:
                chosen = next(
                    (
                        path for path in alternatives
                        if record["id"].split("-")[-1] in path.name.lower()
                    ),
                    alternatives[0],
                )
                contents = chosen.read_text(encoding="latin-1").splitlines()
            if contents is None:
                for archive in sorted(source_directory.glob("*.tar")):
                    with tarfile.open(archive) as package:
                        members = [
                            member for member in package.getmembers()
                            if Path(member.name).name.lower().startswith(name.lower())
                        ]
                        preferred = [
                            member for member in members
                            if (
                                ("noso" in record["id"] and "noso" in member.name)
                                or (
                                    record["id"].endswith("-so")
                                    and member.name.endswith("-so")
                                )
                            )
                        ]
                        selected = (preferred or members)
                        if selected:
                            extracted = package.extractfile(selected[0])
                            if extracted is not None:
                                contents = extracted.read().decode(
                                    "latin-1"
                                ).splitlines()
                                break
            if contents is None:
                for archive in sorted(source_directory.glob("*.zip")):
                    with zipfile.ZipFile(archive) as package:
                        members = [
                            member
                            for member in package.namelist()
                            if Path(member).name.lower().startswith(name.lower())
                        ]
                        if members:
                            contents = package.read(members[0]).decode(
                                "latin-1"
                            ).splitlines()
                            break
        if contents is None:
            continue
        helper = f"potlibc_embedded_data_{number}"

        def replace(match: re.Match[str]) -> str:
            arguments = match.group(2)
            unit = re.search(r"(?i)\bunit\s*=\s*([^,]+)", arguments)
            unit_expression = (
                unit.group(1).strip() if unit else arguments.split(",", 1)[0].strip()
            )
            indent = match.group(1)
            return (
                f"{indent}close(unit={unit_expression},status='delete')\n"
                f"{indent}open(unit={unit_expression},status='scratch',"
                "action='readwrite')\n"
                f"{indent}call {helper}({unit_expression})\n"
                f"{indent}rewind({unit_expression})"
            )

        text = pattern.sub(replace, text)
        writes = "\n".join(
            f"  write(unit,'(a)') '{line.replace(chr(39), chr(39) * 2)}'"
            for line in contents
        )
        helpers.append(
            f"\nsubroutine {helper}(unit)\n"
            "  implicit none\n"
            "  integer, intent(in) :: unit\n"
            f"{writes}\n"
            f"end subroutine {helper}\n"
        )
    return text + "".join(helpers)


def _standalone_program_module(text: str, entry_name: str) -> str:
    program = re.search(r"(?im)^\s*program\s+\w+\s*$", text)
    contains = re.search(r"(?im)^\s*contains\s*$", text)
    if not program or not contains:
        return text
    prefix = text[:program.start()]
    body = text[program.end():contains.start()]
    lines = body.splitlines()
    executable = next(
        (
            index
            for index, line in enumerate(lines)
            if re.match(
                r"(?i)^\s*(?:open|read|write|print|call|do\b|"
                r"[a-z_]\w*(?:\([^=]*\))?\s*=)",
                line,
            )
        ),
        None,
    )
    if executable is None:
        raise ValueError(f"{entry_name}: standalone driver has no initializer")
    declarations = "\n".join(lines[:executable])
    def initialize_real_declaration(match: re.Match[str]) -> str:
        values = re.split(r",(?![^()]*\))", match.group(2))
        return match.group(1) + ",".join(
            f"{value.strip()}=0.0d0" for value in values
        )

    declarations = re.sub(
        r"(?im)^(\s*double\s+precision\s*::\s*)([^\n=]+?)\s*$",
        initialize_real_declaration,
        declarations,
    )
    initialization = "\n".join(
        line for line in lines[executable:]
        if not re.match(r"(?i)^\s*print\s*", line)
    )
    procedures = text[contains.end():]
    procedures = re.sub(
        r"(?im)^\s*end\s*(?:program(?:\s+\w+)?)?\s*$",
        "end module potlibc_standalone_pes",
        procedures,
        count=1,
    )
    return (
        prefix
        + "module potlibc_standalone_pes\n"
        + declarations
        + "\ncontains\n"
        + "subroutine potlibc_standalone_init()\n"
        + initialization
        + "\nend subroutine potlibc_standalone_init\n"
        + procedures
    )


def _save_entry_locals(text: str) -> str:
    """Ensure program units with ENTRY statements retain initialized locals."""
    unit_pattern = re.compile(
        r"(?ims)(^\s*(?:(?:recursive|pure|elemental|impure|module)\s+)*"
        r"(?:(?:double\s+precision|real(?:\s*\*\s*\d+|\s*\([^)\n]*\))?|"
        r"integer(?:\s*\*\s*\d+|\s*\([^)\n]*\))?|"
        r"logical(?:\s*\*\s*\d+|\s*\([^)\n]*\))?|"
        r"character(?:\s*\*\s*(?:\d+|\([^)\n]*\))|\s*\([^)\n]*\))?|"
        r"complex(?:\s*\*\s*\d+|\s*\([^)\n]*\))?|"
        r"type\s*\(\s*dual2\s*\))\s+)?(?:subroutine|function)\b.*?"
        r"^\s*end(?:\s+(?:subroutine|function)(?:\s+\w+)?)?\s*$)"
    )

    def replace(match: re.Match[str]) -> str:
        unit = match.group(1)
        if not re.search(r"(?im)^\s*entry\b", unit):
            return unit
        lines = unit.splitlines()
        if any(re.match(r"(?i)^\s*save\s*$", line) for line in lines):
            return unit
        cleaned = [
            line for line in lines if not re.match(r"(?i)^\s*save\b", line)
        ]
        header = next(
            index
            for index, line in enumerate(cleaned)
            if re.search(r"(?i)\b(?:subroutine|function)\b", line)
        )
        insert_at = header + 1
        for index, line in enumerate(cleaned[insert_at:], insert_at):
            if re.match(r"(?i)^\s*(?:use\s+\w+|implicit\b)", line):
                insert_at = index + 1
            elif line.strip() and not line.lstrip().startswith("!"):
                break
        cleaned.insert(insert_at, "save")
        return "\n".join(cleaned)

    return unit_pattern.sub(replace, text)


def _read_real_components(text: str) -> str:
    """Read numeric records into the real lane of transformed dual variables."""
    nondual: set[str] = set()
    declaration = re.compile(
        r"(?i)^[ \t]*(?:real(?:\s*\([^)\n]*\)|\s*\*\s*\d+)?|"
        r"double[ \t]+precision|integer|logical|character(?:\s*\([^)\n]*\)|"
        r"\s*\*\s*\d+)?)\b(.*)$"
    )
    for line in text.splitlines():
        item = declaration.match(line)
        if item is None:
            continue
        variables = item.group(1)
        if "::" in variables:
            variables = variables.split("::", 1)[1]
        for name in re.findall(r"(?i)\b([a-z]\w*)\b", variables):
            nondual.add(name.lower())
    pattern = re.compile(
        r"(?im)^([ \t]*read[ \t]*\([^)\n]+\))[ \t]*([^ \t\n].*)$"
    )

    def replace(match: re.Match[str]) -> str:
        targets = match.group(2)
        indent = re.match(r"\s*", match.group(1)).group(0)

        implied = re.fullmatch(
            r"(?i)\(\s*([a-z]\w*\s*\([^()]*\))\s*,\s*"
            r"([a-z]\w*)\s*=\s*([^,()]+)\s*,\s*([^,()]+)\s*\)",
            targets,
        )
        if implied is not None:
            target, index, lower, upper = implied.groups()
            name = re.match(r"(?i)\s*([a-z]\w*)", target).group(1)
            if name[0].lower() not in "ijklmn" and name.lower() not in nondual:
                return (
                    f"{indent}block\n"
                    f"{indent}  real*8 potlibc_read_values("
                    f"{lower}:{upper})\n"
                    f"{indent}  {match.group(1).strip()} "
                    f"(potlibc_read_values({index}),{index}={lower},{upper})\n"
                    f"{indent}  do {index}={lower},{upper}\n"
                    f"{indent}    {target}=potlibc_read_values({index})\n"
                    f"{indent}  end do\n"
                    f"{indent}end block"
                )

        assignments: list[tuple[str, str]] = []

        def temporary(token: re.Match[str]) -> str:
            name = token.group(1)
            suffix = token.group(2) or ""
            if name[0].lower() in "ijklmn" or name.lower() in {
                "end",
                "err",
                "iostat",
            } | nondual:
                return token.group(0)
            target = f"{name}{suffix}"
            temporary_name = f"potlibc_read_value_{len(assignments) + 1}"
            assignments.append((target, temporary_name))
            return temporary_name

        converted = re.sub(
            r"(?i)\b([a-z]\w*)(\s*\([^()\n]*\))?(?!\s*%)",
            temporary,
            targets,
        )
        if not assignments:
            return match.group(0)
        declarations = ",".join(item[1] for item in assignments)
        stores = "\n".join(
            f"{indent}  {target}={temporary_name}"
            for target, temporary_name in assignments
        )
        return (
            f"{indent}block\n"
            f"{indent}  real*8 {declarations}\n"
            f"{indent}  {match.group(1).strip()} {converted}\n"
            f"{stores}\n"
            f"{indent}end block"
        )

    return pattern.sub(replace, text)


def _suppress_formatted_output(text: str) -> str:
    character_names: set[str] = set()
    input_units = {
        match.group(1).lower()
        for match in re.finditer(
            r"(?i)\bread\s*\(\s*([a-z]\w*|\d+)", text
        )
    }
    declarations = re.finditer(
        r"(?im)^\s*character(?:\s*\*\s*(?:\(\s*\*\s*\)|\d+)|"
        r"\s*\([^)]*\))?\s*(?:::\s*)?([^!\n]+)",
        text,
    )
    for declaration in declarations:
        for item in declaration.group(1).split(","):
            name = re.match(r"\s*([a-z]\w*)", item, re.I)
            if name:
                character_names.add(name.group(1).lower())

    statement = re.compile(
        r"(?im)^([ \t]*(?:\d+[ \t]+)?)(?:"
        r"write\s*\(\s*([a-z]\w*|\d+|\*)[^\n]*|print\s+[^\n]*)$"
    )

    def replace(match: re.Match[str]) -> str:
        target = match.group(2)
        if target is not None and (
            target.lower() in character_names or target.lower() in input_units
        ):
            return match.group(0)
        return match.group(1) + "continue"

    return statement.sub(replace, text)


def _callable(record: dict[str, Any]) -> tuple[str, tuple[str, ...]]:
    declarations = [
        item
        for item in record["declarations"]
        if item["kind"] in {"subroutine", "function", "entry"}
        and not item["signature"].lower().startswith("end ")
    ]
    by_name = {item["name"].lower(): item for item in declarations}
    preferred = {
        "brh2so": "pot11",
        "bjkkl": "trbunk",
        "oh2dim": "pot2",
        "c2h7pot5": "c2h7pot",
    }.get(record["id"])
    if record["path"].endswith("/no2-pipnn.f"):
        preferred = "evvdvdx"
    item = by_name.get(preferred) if preferred else record.get("entry")
    if item is None:
        raise ValueError(f"{record['id']}: inventory has no callable entry")
    return (
        str(item["name"]).lower(),
        tuple(str(arg).lower().replace(" ", "") for arg in item["arguments"]),
    )


def _plan(record: dict[str, Any], vendor: Path) -> SurfacePlan:
    atom_count = record["atom_count"]
    variable_atoms = atom_count["kind"] == "variable"
    source_directory = (vendor / record["path"]).parent
    package_members: set[str] = set()
    for archive in source_directory.glob("*.tar"):
        with tarfile.open(archive) as package:
            package_members.update(
                Path(member.name).name.lower() for member in package.getmembers()
            )
    for archive in source_directory.glob("*.zip"):
        with zipfile.ZipFile(archive) as package:
            package_members.update(
                Path(member).name.lower() for member in package.namelist()
            )
    missing_parameters = [
        reference["name"]
        for reference in record["external_data_files"]
        if reference.get("literal")
        and not reference["resolved_paths"]
        and not any(source_directory.glob(reference["name"] + ".*"))
        and not any(
            member.startswith(reference["name"].lower())
            for member in package_members
        )
    ]
    if missing_parameters:
        raise ValueError(
            f"{record['id']}: required parameter data is absent from the vendor "
            f"archive: {', '.join(missing_parameters)}"
        )
    atoms = None if variable_atoms else int(atom_count["value"])
    if atoms is not None and atoms < 1:
        raise ValueError(f"{record['id']}: invalid fixed atom count {atoms}")

    entry_name, arguments = _callable(record)
    family = str(record["interface_family"]).upper()
    source = vendor / record["path"]
    text = _implementation_text(record, vendor)
    module_name = _entry_module(text, entry_name)
    common_name = None

    commons = {str(name).lower() for name in record["common_blocks"]}
    if (
        family == POTLIB_FAMILY
        and entry_name == "pot"
        and not arguments
        and "usrocm" in commons
    ):
        strategy = "potlib-common"
    elif (
        (
            family in PAIR_FAMILIES | {"4-XS", "5-1S", "6-1S", "7-1S"}
            or family == POTLIB_FAMILY
        )
        and entry_name == "pot"
        and not arguments
    ):
        if "pt1cm" in commons:
            common_name = "pt1cm"
        elif "potcm" in commons:
            common_name = "potcm"
        else:
            raise ValueError(
                f"{record['id']}: scalar {family} has neither PT1CM nor POTCM"
            )
        strategy = "pair-common"
    elif family in DISTANCE_FAMILIES and entry_name == "pot":
        if not (3 <= len(arguments) <= 5):
            raise ValueError(
                f"{record['id']}: unresolved {family} signature "
                f"{entry_name}{arguments!r}"
            )
        strategy = "distance-batch"
    elif entry_name == "surf" and len(arguments) == 4:
        strategy = "section2-surf"
    elif entry_name == "pot" and _is_cartesian_signature(arguments):
        strategy = "cartesian-pot"
    elif family.startswith("HO-MM-"):
        strategy = "ho-mm"
    elif entry_name in {"c2h7pot", "c2h7poth", "f"}:
        strategy = "flat-cartesian"
    elif entry_name in {"pesh3cl", "ho3_pes_interface", "h3o2pipnn"}:
        strategy = "matrix-cartesian"
    elif entry_name == "fh2onn":
        strategy = "fh2o-nn"
    elif entry_name == "fh2o_soc_nn":
        strategy = "fh2o-soc"
    elif entry_name in {"n3pes", "o3pes"}:
        strategy = "triatomic-function"
    elif entry_name == "o4pes":
        strategy = "distance-function"
    elif entry_name == "n4pipnn":
        strategy = "n4-cartesian"
    elif entry_name == "vibpot":
        strategy = "distance-vector"
    elif entry_name == "energy_gradient":
        strategy = "energy-gradient"
    elif entry_name == "evvdvdx":
        strategy = "no2-nn"
    elif entry_name == "pot11":
        strategy = "distance-pot11"
    elif entry_name == "trbunk":
        strategy = "four-body-trbunk"
    elif entry_name == "pot2":
        strategy = "pot2-common"
    elif entry_name == "pot" and len(arguments) == 9:
        strategy = "multistate-cartesian"
    elif entry_name == "pot" and len(arguments) == 11:
        strategy = "nh3-cartesian"
    elif family == "4-1S" and entry_name == "pot":
        strategy = "four-body-diatom"
    elif entry_name == "pot" and not arguments:
        strategy = "pair-common"
        common_name = "pt1cm"
    else:
        raise ValueError(
            f"{record['id']}: unresolved {family} signature "
            f"{entry_name}{arguments!r}"
        )
    return SurfacePlan(
        record, source, atoms, variable_atoms, strategy, entry_name, arguments,
        module_name, common_name
    )


def _is_cartesian_signature(arguments: tuple[str, ...]) -> bool:
    if len(arguments) != 7:
        return False
    compact = tuple(re.sub(r"_$", "", arg) for arg in arguments)
    return (
        compact[:4] == ("x", "y", "z", "e")
        and compact[4:] in {
            ("dedx", "dedy", "dedz"),
            ("dvdx", "dvdy", "dvdz"),
        }
    )


def _transform(plan: SurfacePlan, vendor: Path) -> tuple[str, str]:
    original = _implementation_text(plan.record, vendor)
    if plan.record["id"] == "n4pes-l-imls-g2-v1":
        package_path = plan.source.with_suffix(".tar")
        prepot_member = f"{package_path.stem}/prepot.f"
        with tarfile.open(package_path) as package:
            extracted = package.extractfile(prepot_member)
            if extracted is None:
                raise ValueError(
                    f"{plan.record['id']}: archive has no {prepot_member}"
                )
            prepot = extracted.read().decode("latin-1")
        prepot = re.sub(
            r"(?im)^(\s*subroutine\s+readpes\s*)$",
            r"\1\n      use potlibc_n4_data_state",
            prepot,
        )
        prepot = re.sub(
            r"(?i)(\b(?:file|inquire\s*\(\s*file)\s*=\s*)['\"]pes\.dat['\"]",
            r"\1potlibc_n4_path",
            prepot,
        )
        original = """\
      module potlibc_n4_data_state
      use iso_c_binding, only: c_int
      implicit none
      character(128) :: potlibc_n4_path
      interface
      integer(c_int) function potlibc_n4_data_fd() bind(c)
      import c_int
      end function
      end interface
      contains
      subroutine potlibc_n4_prepare()
      integer(c_int) :: fd
      fd=potlibc_n4_data_fd()
      write(potlibc_n4_path,'("/proc/self/fd/",i0)') fd
      end subroutine
      end module potlibc_n4_data_state
""" + original + "\n" + prepot
    original = re.sub(
        r"(?im)^(\s*read\s*\([^)\n]+\))\s*comlin\s*$",
        r"\1",
        original,
    )
    if plan.record["id"] == "bjkkl":
        original = re.sub(
            r"(?ims)^\s*program\s+dummy\b.*?^\s*end\s*$",
            "",
            original,
            count=1,
        )
        original = re.sub(
            r"(?i)zthreej\(([^)\n]*),\s*0(?:\.0*)?d0\s*\)",
            r"zthreej(\1,dual2(0d0,0d0,0d0,0d0))",
            original,
        )
    if plan.strategy in {"triatomic-function", "distance-function"}:
        original = _standalone_program_module(original, plan.entry_name)
    fixed_unit_init = (
        _active_call_count(original, "potinfo") == 1
        and _active_call_count(original, "ancvrt") == 1
    )
    cbe = plan.record["path"] == "src/CH5/ch52008.f"
    real_source = adapted_source(original, fixed_unit_init, cbe)
    real_source = re.sub(
        r"(?im)^(\s*)real\s+function\b",
        r"\1double precision function",
        real_source,
    )
    real_source = re.sub(r"(?i)\bidnint\b", "nint", real_source)
    real_source = re.sub(r"(?i)\bdatan2d\b", "atan2", real_source)
    real_source = re.sub(
        r"(?i)\bisnan\s*\(\s*([a-z_][a-z0-9_]*)\s*\)",
        r"(\1.ne.\1)",
        real_source,
    )
    real_source = re.sub(
        r"(?im)^(\s*)deep\s+neural\s+network\s*$",
        r"\1! deep neural network",
        real_source,
    )
    real_source = re.sub(
        r"(?i)maxval\s*\(\s*([^,()]+)\s*,\s*nbonds\s*\)",
        r"maxval(\1)",
        real_source,
    )
    real_source = re.sub(
        r"(?i)minval\s*\(\s*rb\s*(?:\(\s*:\s*\))?\s*\)",
        "min(min(rb(1),rb(2)),min(rb(3),rb(4)))",
        real_source,
    )
    real_source = re.sub(
        r"(?i)sum\s*\(\s*(vec\s*\(\s*:\s*,\s*\d+\s*\)|"
        r"crop[12]\s*\(\s*:\s*\))\s*\*\*\s*2\s*\)",
        lambda match: "+".join(
            f"{match.group(1).replace(':', str(index))}**2"
            for index in range(1, 4)
        ),
        real_source,
    )
    for include in re.findall(r"(?im)^\s*include\s*['\"]([^'\"]+)['\"]", real_source):
        candidates = sorted(plan.source.parent.rglob(Path(include).name))
        if not candidates:
            contents = ""
            for archive in sorted(plan.source.parent.glob("*.tar")):
                with tarfile.open(archive) as package:
                    member = next(
                        (
                            item for item in package.getmembers()
                            if Path(item.name).name == Path(include).name
                        ),
                        None,
                    )
                    extracted = package.extractfile(member) if member else None
                    if extracted is not None:
                        contents = extracted.read().decode("latin-1")
                        break
            if not contents:
                first_common = re.search(r"(?im)^\s*common\s*/", real_source)
                first_executable = (
                    re.search(
                        r"(?im)^\s*(?:do|if|call|open|read|[a-z]\w*\s*=)\b",
                        real_source[first_common.start():],
                    )
                    if first_common
                    else None
                )
                if first_common and first_executable:
                    declaration_end = first_common.start() + first_executable.start()
                    contents = real_source[first_common.start():declaration_end]
                else:
                    raise ValueError(
                        f"{plan.record['id']}: include {include!r} is absent and "
                        "no declaration preamble can replace it"
                    )
        elif len(candidates) == 1:
            contents = _source_text(candidates[0])
        else:
            raise ValueError(
                f"{plan.record['id']}: include {include!r} resolves to "
                f"{len(candidates)} files"
            )
        real_source = re.sub(
            rf"(?im)^\s*include\s*['\"]{re.escape(include)}['\"]\s*$",
            contents,
            real_source,
        )
    dependencies: list[str] = []
    if plan.record["id"] in {"pes-shell", "soc-fh2o-nn"}:
        for dependency in plan.record["source_dependencies"]:
            dependency_path = vendor / dependency["path"]
            dependencies.append(
                _free_source(_source_text(dependency_path), dependency_path.suffix)
            )
    real_source = "\n".join(dependencies) + _free_source(
        real_source, plan.source.suffix, normalize_source=False
    )
    real_source = _suppress_formatted_output(real_source)
    real_source = re.sub(
        r"(?i)minval\s*\(\s*rb\s*(?:\(\s*:\s*\))?\s*\)",
        "min(min(rb(1),rb(2)),min(rb(3),rb(4)))",
        real_source,
    )
    real_source = re.sub(
        r"(?i)sum\s*\(\s*(vec\s*\(\s*:\s*,\s*\d+\s*\)|"
        r"crop[12]\s*\(\s*:\s*\))\s*\*\*\s*2\s*\)",
        lambda match: "+".join(
            f"{match.group(1).replace(':', str(index))}**2"
            for index in range(1, 4)
        ),
        real_source,
    )
    if plan.strategy in {
        "flat-cartesian", "matrix-cartesian", "fh2o-nn", "fh2o-soc"
    }:
        real_source = re.sub(
            r"(?im)^(\s*)real(?=\s*(?:::|,))", r"\1REAL*8", real_source
        )
    real_source = re.sub(r"(?i)\breal\*8\s*::", "REAL*8 ", real_source)
    real_source = re.sub(
        r"(?i)min\s*\(([^,()\n]+),([^,()\n]+),([^,()\n]+)\)",
        r"min(min(\1,\2),\3)",
        real_source,
    )
    real_source = re.sub(
        r"(?im)^(\s*data\s+tiny\s*/\s*)([-+]?(?:\d+\.\d*|\.\d+)"
        r"(?:[dD][+-]?\d+)?)(\s*/)",
        r"\1dual2(\2,0d0,0d0,0d0)\3",
        real_source,
    )
    if plan.record["id"] in {"nich4", "niheam5", "niheam6", "niheam7"}:
        archive = plan.record["archive_source"]
        package_path = vendor / archive["path"]
        with tarfile.open(package_path) as package:
            archive_directory = Path(archive["member"]).parent
            parameter_members = sorted(
                [
                member
                for member in package.getmembers()
                if Path(member.name).parent == archive_directory
                and re.fullmatch(r"fort\.(\d+)", Path(member.name).name)
                ],
                key=lambda member: member.name,
            )
            helpers: list[str] = []
            for member in parameter_members:
                extracted = package.extractfile(member)
                if extracted is None:
                    raise ValueError(
                        f"{plan.record['id']}: cannot extract {member.name}"
                    )
                parameter_lines = (
                    extracted.read().decode("latin-1").splitlines()
                )
                writes = "\n".join(
                    f"  write(unit,'(a)') "
                    f"'{line.replace(chr(39), chr(39) * 2)}'"
                    for line in parameter_lines
                )
                unit_number = Path(member.name).suffix[1:]
                helpers.append(
                    f"\nsubroutine potlibc_archived_unit{unit_number}(unit)\n"
                    "  implicit none\n"
                    "  integer, intent(in) :: unit\n"
                    f"{writes}\n"
                    f"end subroutine potlibc_archived_unit{unit_number}\n"
                )
        real_source = re.sub(
            r"(?im)^(\s*)read\s*\(\s*4\s*,",
            r"\1rewind(4)\n\1read(4,",
            real_source,
            count=1,
        )
        real_source += "".join(helpers)
    real_source = _embed_parameter_files(real_source, plan.record, vendor)
    if plan.record["id"] == "soc-fh2o-nn":
        dependency_directory = vendor / "src/done/FH2O"
        dependency_record = {
            "id": plan.record["id"],
            "path": "src/done/FH2O/fh2o-pipnn.f",
            "external_data_files": [
                {
                    "name": name,
                    "literal": True,
                    "resolved_paths": [
                        (dependency_directory / name).relative_to(vendor).as_posix()
                    ],
                }
                for name in ("biases.txt", "weights.txt")
            ],
        }
        real_source = _embed_parameter_files(
            real_source, dependency_record, vendor
        )

    if plan.strategy in {"potlib-common", "pair-common", "pot2-common"}:
        utility_path = vendor / "utilities" / "utility.f"
        real_utility = _free_source(_source_text(utility_path), utility_path.suffix)
    else:
        real_utility = (
            "subroutine potlibc_generated_empty_utility()\n"
            "  implicit none\n"
            "end subroutine potlibc_generated_empty_utility\n"
        )

    source_ids, source_commons = harvest(real_source)
    utility_ids, utility_commons = harvest(real_utility)
    source_ids |= set(
        re.findall(
            r"(?im)^\s*entry\s+([a-z_][a-z0-9_]*)", real_source
        )
    )
    identifiers = source_ids | utility_ids
    external_names = {
        "to", "dsyev", "dgemm", "dgemv", "dgesv", "dgetrf", "dgetri",
        "exit", "stop",
    }
    identifiers = {
        identifier for identifier in identifiers
        if identifier.lower() not in external_names
    }
    commons = source_commons | utility_commons
    transformed_source = to_dual2(real_source, identifiers, commons)
    transformed_source = re.sub(
        r"(?i)\bdsyev\b", "potlibc_dsyev", transformed_source
    )
    transformed_source = _read_real_components(transformed_source)
    transformed_source = _save_entry_locals(transformed_source)
    if plan.record["id"] == "bjkkl":
        transformed_source = re.sub(
            r"(?i)zthreej\(([^)\n]*),\s*0(?:\.0*)?d0\s*\)",
            r"zthreej(\1,dual2(0d0,0d0,0d0,0d0))",
            transformed_source,
        )
    if plan.record["id"] == "h3-bh":
        transformed_source = re.sub(
            r"(?i)call\s+zgolden\s*\(\s*1\.3d0\s*,\s*1\.4d0\s*,"
            r"\s*1\.5d0\s*,\s*1\.d-14\s*,",
            "call zgolden(dual2(1.3d0,0d0,0d0,0d0),"
            "dual2(1.4d0,0d0,0d0,0d0),"
            "dual2(1.5d0,0d0,0d0,0d0),"
            "dual2(1.d-14,0d0,0d0,0d0),",
            transformed_source,
        )
    transformed_source = re.sub(
        r"dual2\(dual2\(([^,()\n]+),0d0,0d0,0d0\),0d0,0d0,0d0\)",
        r"dual2(\1,0d0,0d0,0d0)",
        transformed_source,
        flags=re.I,
    )
    if (
        plan.module_name
        and re.search(
        r"(?im)^\s*subroutine\s+zprepot\b", transformed_source
        )
        and not re.search(
            r"(?im)^\s*public\s*::[^\n]*\bzprepot\b", transformed_source
        )
    ):
        transformed_source, count = re.subn(
            r"(?im)^(\s*public\s*::[^\n]*\bzpot\b[^\n]*)$",
            r"\1, zprepot",
            transformed_source,
            count=1,
        )
        if count != 1:
            raise ValueError(
                f"{plan.record['id']}: module PREPOT cannot be made public"
            )
    transformed_source = re.sub(
        r"(,\s*zprepot)(?:,\s*zprepot)+",
        r"\1",
        transformed_source,
        flags=re.I,
    )
    transformed_utility = to_dual2(real_utility, identifiers, commons)
    transformed_utility = re.sub(
        r"(?i)\bdsyev\b", "potlibc_dsyev", transformed_utility
    )
    transformed_utility = _read_real_components(transformed_utility)
    if cbe:
        transformed_source = real_source + "\n" + transformed_source
        transformed_utility = real_utility + "\n" + transformed_utility
    return transformed_source, transformed_utility


def _distance_call(plan: SurfacePlan) -> str:
    arguments = plan.entry_arguments
    triatomic_order = _internal_distance_order(plan)
    if triatomic_order is None:
        raise ValueError(
            f"{plan.record['id']}: cannot derive triatomic coordinate order"
        )
    pair_value = {
        (1, 2): "distance3(coords(1,:),coords(2,:))",
        (2, 3): "distance3(coords(2,:),coords(3,:))",
        (1, 3): "distance3(coords(1,:),coords(3,:))",
    }
    declarations = [
        "  type(dual2) :: r(1,3), energies(1), derivatives(3,1)",
        *(
            f"  r(1,{index})={pair_value[pair]}"
            for index, pair in enumerate(triatomic_order, 1)
        ),
    ]
    actual: list[str] = []
    for index, argument in enumerate(arguments):
        if index == 0:
            actual.append("r")
        elif index == 1:
            actual.append("energies")
        elif argument in {"de", "der", "dedr", "derivatives"}:
            actual.append("derivatives")
        elif argument in {"nt", "n", "ni", "nvals"}:
            actual.append("1")
        elif argument in {"nsurf", "ndum"}:
            actual.append("1")
        else:
            raise ValueError(
                f"{plan.record['id']}: unresolved distance argument {argument!r}"
            )
    declarations += [f"  call zpot({','.join(actual)})", "  value=energies(1)"]
    return "\n".join(declarations)


def _initialization_name(plan: SurfacePlan) -> str | None:
    if plan.record["id"] == "n4pes-l-imls-g2-v1":
        return "prepot"
    if plan.record["id"] in {"h3pluskbnnd", "h3plusunod"}:
        return None
    declarations = {
        item["name"].lower()
        for item in plan.record["declarations"]
        if item["kind"] in {"subroutine", "entry"}
    }
    if plan.strategy in {"triatomic-function", "distance-function"}:
        return "potlibc_standalone_init"
    name = next(
        (
            name
            for name in ("prepot", "pes_init", "bjinit", "prepot11")
            if name in declarations
        ),
        None,
    )
    if name == "prepot":
        source = _source_text(plan.source)
        entry = re.search(r"(?im)^[ \t]*entry\b", source)
        if entry is not None and not re.search(
            r"(?im)^[ \t]*(?:\d+[ \t]+)?return\b", source[: entry.start()]
        ):
            return None
    return name


def _normalize_label(raw: str) -> str:
    label = raw.strip().lower()
    label = label.replace("first", "1").replace("second", "2").replace(
        "third", "3"
    )
    label = label.replace("'", "p")
    label = re.sub(r"[^a-z0-9]+", "", label)
    return label


def _metadata_symbols(plan: SurfacePlan) -> list[str] | None:
    metadata = plan.record.get("atomic_metadata", {})
    if metadata.get("kind") != "fixed":
        return None
    symbols = metadata.get("symbols")
    if isinstance(symbols, list) and len(symbols) == 3:
        return [str(symbol).lower() for symbol in symbols]
    numbers = metadata.get("numbers")
    if isinstance(numbers, list) and len(numbers) == 3:
        return [f"z{int(number)}" for number in numbers]
    return None


def _label_symbol(label: str, symbols: list[str]) -> str | None:
    unqualified = re.sub(r"^[123]", "", label)
    unqualified = re.sub(r"p$", "", unqualified)
    matches = [
        symbol for symbol in set(symbols)
        if unqualified == symbol
    ]
    return matches[0] if len(matches) == 1 else None


def _label_rank(label: str) -> tuple[int, str]:
    if label.startswith(("1",)):
        return (0, label)
    if label.startswith(("2",)) or label.endswith("p"):
        return (1, label)
    if label.startswith(("3",)):
        return (2, label)
    return (0, label)


def _internal_distance_order(plan: SurfacePlan) -> tuple[tuple[int, int], ...] | None:
    """Parse header R(1..3) definitions for triatomic internal coordinates."""
    atom_count = plan.record.get("atom_count", {})
    if atom_count.get("kind") != "fixed" or int(atom_count.get("value", 0)) != 3:
        return None
    symbols = _metadata_symbols(plan)
    source = _source_text(plan.source)
    coordinates = str(plan.record.get("documentation", {}).get("coordinates", ""))
    text = source[:14000] + "\n" + coordinates
    definitions: dict[int, tuple[str, str]] = {}
    for match in re.finditer(
        r"(?i)\br\s*\(\s*([123])\s*\)\s*=\s*r\s*\(\s*([^)]+)\s*\)",
        text,
    ):
        index = int(match.group(1))
        if index in definitions:
            continue
        descriptor = match.group(2)
        if "-" not in descriptor:
            continue
        left, right = (part.strip() for part in descriptor.split("-", 1))
        definitions[index] = (_normalize_label(left), _normalize_label(right))
    if len(definitions) != 3:
        if symbols is not None:
            raise ValueError(
                f"{plan.record['id']}: fixed triatomic metadata requires "
                "three internal-coordinate bond definitions"
            )
        if plan.strategy in {
            "pair-common",
            "pot2-common",
            "distance-batch",
            "distance-pot11",
        }:
            return ((1, 2), (2, 3), (1, 3))
        return None

    by_index = [definitions[1], definitions[2], definitions[3]]
    if symbols is not None:
        labels = {label for pair in by_index for label in pair if label}
        positions: dict[str, list[int]] = {}
        for index, symbol in enumerate(symbols, 1):
            positions.setdefault(symbol, []).append(index)
        label_groups: dict[str, list[str]] = {}
        for label in labels:
            symbol = _label_symbol(label, symbols)
            if symbol is None:
                break
            label_groups.setdefault(symbol, []).append(label)
        else:
            if all(
                len(label_groups.get(symbol, [])) == len(indexes)
                for symbol, indexes in positions.items()
            ):
                order: dict[str, int] = {}
                for symbol, indexes in positions.items():
                    ordered_labels = sorted(label_groups[symbol], key=_label_rank)
                    order.update(zip(ordered_labels, indexes))
                return tuple(
                    tuple(sorted((order[left], order[right])))
                    for left, right in by_index
                )

    sets = [set(pair) for pair in by_index]
    tokens = sorted({label for pair in by_index for label in pair if label})
    if len(tokens) == 3 and all(len(item) == 2 for item in sets):
        first_atom = sets[0] & sets[2]
        second_atom = sets[0] & sets[1]
        third_atom = sets[1] & sets[2]
        if (
            len(first_atom) == 1
            and len(second_atom) == 1
            and len(third_atom) == 1
        ):
            order = {
                next(iter(first_atom)): 1,
                next(iter(second_atom)): 2,
                next(iter(third_atom)): 3,
            }
            return tuple(
                tuple(sorted((order[left], order[right])))
                for left, right in by_index
            )
    if len(tokens) >= 3:
        order = {token: index + 1 for index, token in enumerate(tokens[:3])}
        pairs: list[tuple[int, int]] = []
        for left, right in by_index:
            if left not in order or right not in order or left == right:
                return None
            pair = tuple(sorted((order[left], order[right])))
            if pair not in {(1, 2), (2, 3), (1, 3)}:
                return None
            pairs.append(pair)
        if len(set(pairs)) == 3:
            return tuple(pairs)
        return None

    if len(tokens) == 2:
        repeated = next((pair for pair in by_index if pair[0] == pair[1]), None)
        if repeated is None:
            return None
        mixed_indexes = [index for index, pair in enumerate(by_index) if pair[0] != pair[1]]
        if len(mixed_indexes) != 2:
            return None
        if symbols is not None and symbols[0] == symbols[2]:
            repeated_pair = (1, 3)
            mapping_by_position = {
                mixed_indexes[0]: (1, 2),
                mixed_indexes[1]: (2, 3),
            }
        else:
            repeated_pair = (2, 3)
            mapping_by_position = {
                mixed_indexes[0]: (1, 2),
                mixed_indexes[1]: (1, 3),
            }
        pairs: list[tuple[int, int]] = []
        for index, pair in enumerate(by_index):
            if pair[0] == pair[1]:
                pairs.append(repeated_pair)
            else:
                pairs.append(mapping_by_position[index])
        return tuple(pairs)
    return None


# CODATA 2018. The c2h6 adapters already use the Bohr factor below.
BOHR_TO_ANGSTROM = "0.529177210903d0"
KCAL_PER_HARTREE = "627.5094740631d0"
EV_PER_HARTREE = "27.211386245988d0"


def _energy_bridge(plan: SurfacePlan) -> str:
    entry = "z" + plan.entry_name
    if plan.strategy in {"triatomic-function", "distance-function"}:
        module_use = f"  use potlibc_standalone_pes, only: {entry}\n"
    else:
        module_use = (
            f"  use {plan.module_name}, only: {entry}\n"
            if plan.module_name
            else ""
        )
    common = plan.common_name or "pt1cm"
    if plan.record["id"].startswith("c2h6"):
        common_coordinates = 3 * int(plan.atoms)
        body = f"""\
  integer, parameter :: maxatom=25,isurf=5,jsurf=15
  type(dual2) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual2) :: dgscart(maxatom,3),descart(maxatom,3,isurf)
  type(dual2) :: dijcart(maxatom,3,jsurf),cart(maxatom,3),anuzero
  type(dual2) :: q({common_coordinates}),pdot({common_coordinates})
  integer :: nullbl(maxatom),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: atom,axis,k
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  common /zqpdot/ q,pdot
  do atom=1,natoms
    do axis=1,3
      k=3*(atom-1)+axis
      q(k)=coords(atom,axis)*{BOHR_TO_ANGSTROM}
    end do
  end do
  nder=0
  call zpot()
  value=pengygs"""
    elif plan.strategy == "potlib-common":
        body = """\
  integer, parameter :: maxatom=25,isurf=5,jsurf=15
  type(dual2) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual2) :: dgscart(maxatom,3),descart(maxatom,3,isurf)
  type(dual2) :: dijcart(maxatom,3,jsurf),cart(maxatom,3),anuzero
  integer :: nullbl(maxatom),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: atom,axis
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do atom=1,natoms
    do axis=1,3
      cart(atom,axis)=coords(atom,axis)
    end do
  end do
  call zpot()
  value=pengygs"""
    elif plan.strategy in {"pair-common", "pot2-common"}:
        pair_count = int(plan.atoms) * (int(plan.atoms) - 1) // 2
        triatomic_order = _internal_distance_order(plan)
        if pair_count == 3 and triatomic_order is not None:
            pair_value = {
                (1, 2): "distance3(coords(1,:),coords(2,:))",
                (2, 3): "distance3(coords(2,:),coords(3,:))",
                (1, 3): "distance3(coords(1,:),coords(3,:))",
            }
            assignments = "\n".join(
                f"  r({index})={pair_value[pair]}"
                for index, pair in enumerate(triatomic_order, 1)
            )
            body = f"""\
  type(dual2) :: r({pair_count}),energy
  type(dual2) :: dedr({pair_count})
  common /z{common}/ r,energy,dedr
{assignments}
  call {entry}()
  value=energy"""
        else:
            body = f"""\
  type(dual2) :: r({pair_count}),energy
  type(dual2) :: dedr({pair_count})
  integer :: left,right,k
  common /z{common}/ r,energy,dedr
  k=0
  do left=1,natoms-1
    do right=left+1,natoms
      k=k+1
      r(k)=distance3(coords(left,:),coords(right,:))
    end do
  end do
  call {entry}()
  value=energy"""
    elif plan.strategy == "distance-batch":
        body = _distance_call(plan)
    elif plan.strategy == "section2-surf":
        body = f"""\
  type(dual2) :: flat(3*natoms),gradient(3*natoms)
  integer :: atom,axis,k
  do atom=1,natoms
    do axis=1,3
      k=3*(atom-1)+axis
      flat(k)=coords(atom,axis)
    end do
  end do
  call {entry}(value,flat,gradient,3*natoms)"""
    elif plan.strategy == "cartesian-pot":
        body = """\
  type(dual2) :: x(natoms),y(natoms),z(natoms)
  type(dual2) :: dx(natoms),dy(natoms),dz(natoms)
  x=coords(:,1); y=coords(:,2); z=coords(:,3)
  call zpot(x,y,z,value,dx,dy,dz)"""
    elif plan.strategy == "ho-mm":
        if plan.entry_name == "siohpot":
            body = """\
  type(dual2) :: x(natoms),y(natoms),z(natoms),q(natoms),v
  type(dual2) :: dx(natoms),dy(natoms),dz(natoms)
  integer :: indmm(natoms),atom
  x=coords(:,1); y=coords(:,2); z=coords(:,3)
  do atom=1,natoms
    select case (mod(atom-1,3))
    case (0)
      indmm(atom)=14; q(atom)=4
    case (1)
      indmm(atom)=8; q(atom)=-2
    case (2)
      indmm(atom)=1; q(atom)=1
    end select
  end do
  call zsiohpot(x,y,z,q,indmm,v,dx,dy,dz,natoms,natoms,0)
  value=v"""
        elif len(plan.entry_arguments) == 7:
            body = f"""\
  type(dual2) :: x(natoms),y(natoms),z(natoms)
  character(2) :: symbol(natoms)
  x=coords(:,1); y=coords(:,2); z=coords(:,3)
  symbol='AL'
  call {entry}(x,y,z,value,natoms,natoms,symbol)"""
        elif len(plan.entry_arguments) == 9:
            body = f"""\
  type(dual2) :: x(natoms),y(natoms),z(natoms)
  type(dual2) :: dx(natoms),dy(natoms),dz(natoms)
  x=coords(:,1); y=coords(:,2); z=coords(:,3)
  call {entry}(x,y,z,value,dx,dy,dz,natoms,natoms)"""
        else:
            body = f"""\
  type(dual2) :: x(natoms),y(natoms),z(natoms)
  x=coords(:,1); y=coords(:,2); z=coords(:,3)
  call {entry}(x,y,z,value,natoms,natoms)"""
    elif plan.strategy == "flat-cartesian":
        if plan.entry_name == "f":
            body = f"""\
  type(dual2) :: matrix(3,natoms)
  matrix=transpose(coords)*{BOHR_TO_ANGSTROM}
  value=zf(matrix)"""
        elif plan.entry_name == "c2h7pot":
            body = f"""\
  type(dual2) :: xyz(natoms,3)
  xyz=coords*{BOHR_TO_ANGSTROM}
  call {entry}(xyz,value)
  value=value/{KCAL_PER_HARTREE}"""
        else:
            raise AssertionError(plan.entry_name)
    elif plan.strategy == "matrix-cartesian":
        body = f"""\
  type(dual2) :: matrix(3,natoms)
  matrix=transpose(coords)*{BOHR_TO_ANGSTROM}
  call {entry}(matrix,value)
  value=value/{EV_PER_HARTREE}"""
    elif plan.strategy == "fh2o-nn":
        body = f"""\
  type(dual2) :: matrix(3,natoms),states(3)
  matrix=transpose(coords)*{BOHR_TO_ANGSTROM}
  call zfh2onn(matrix,value,states(1),states(2),states(3))
  value=value/{EV_PER_HARTREE}"""
    elif plan.strategy == "fh2o-soc":
        body = f"""\
  type(dual2) :: matrix(3,natoms),vg,socg,vgsoc
  matrix=transpose(coords)*{BOHR_TO_ANGSTROM}
  call zfh2o_soc_nn(matrix,vg,socg,vgsoc)
  value=vg/{KCAL_PER_HARTREE}"""
    elif plan.strategy == "triatomic-function":
        body = f"""\
  type(dual2) :: r12,r13,r23
  r12=distance3(coords(1,:),coords(2,:))*{BOHR_TO_ANGSTROM}
  r13=distance3(coords(1,:),coords(3,:))*{BOHR_TO_ANGSTROM}
  r23=distance3(coords(2,:),coords(3,:))*{BOHR_TO_ANGSTROM}
  value={entry}(r12,r13,r23)/{KCAL_PER_HARTREE}"""
    elif plan.strategy == "distance-function":
        body = f"""\
  type(dual2) :: r(6)
  call pair_distances(coords,natoms,r)
  r=r*{BOHR_TO_ANGSTROM}
  value={entry}(r)/{KCAL_PER_HARTREE}"""
    elif plan.strategy == "distance-vector":
        body = f"""\
  type(dual2) :: r(natoms*(natoms-1)/2)
  call pair_distances(coords,natoms,r)
  call {entry}(r,value,size(r))"""
    elif plan.strategy == "n4-cartesian":
        body = """\
  type(dual2) :: x(natoms),y(natoms),z(natoms)
  type(dual2) :: dx(natoms),dy(natoms),dz(natoms)
  x=coords(:,1); y=coords(:,2); z=coords(:,3)
  call zn4pipnn(x,y,z,value,dx,dy,dz,0)"""
    elif plan.strategy == "energy-gradient":
        states = int(plan.record["id"].split("-")[1])
        body = f"""\
  type(dual2) :: energies({states}),gradient(3,natoms,{states})
  type(dual2) :: coupling({states},{states},natoms,3)
  type(dual2) :: transform({states},{states})
  call zenergy_gradient(coords,energies,gradient,coupling,transform,natoms,{states})
  value=energies(1)"""
    elif plan.strategy == "no2-nn":
        body = f"""\
  type(dual2) :: xcart(3,natoms),vpes(3),dvdxa(3*natoms)
  type(dual2) :: dvdx(3*natoms,3)
  xcart=transpose(coords)*{BOHR_TO_ANGSTROM}
  call zevvdvdx(xcart,value,vpes,dvdxa,dvdx,0)
  value=value/{EV_PER_HARTREE}"""
    elif plan.strategy == "distance-pot11":
        body = """\
  type(dual2) :: r(1,3),energies(1)
  r(1,1)=distance3(coords(1,:),coords(2,:))
  r(1,2)=distance3(coords(2,:),coords(3,:))
  r(1,3)=distance3(coords(1,:),coords(3,:))
  call zpot11(r,energies,1)
  value=energies(1)"""
    elif plan.strategy == "four-body-trbunk":
        body = """\
  type(dual2) :: r12,r34,centers,energy
  r12=distance3(coords(1,:),coords(2,:))
  r34=distance3(coords(3,:),coords(4,:))
  centers=distance3((coords(1,:)+coords(2,:))/2, &
                    (coords(3,:)+coords(4,:))/2)
  call ztrbunk(centers,r12,r34,energy,0,0,0)
  value=energy"""
    elif plan.strategy == "four-body-diatom":
        body = """\
  type(dual2) :: r1,r2,r,th1,th2,tau
  type(dual2) :: a(3),b(3),c(3)
  a=coords(2,:)-coords(1,:); b=coords(4,:)-coords(3,:)
  c=(coords(3,:)+coords(4,:)-coords(1,:)-coords(2,:))/2
  r1=norm3(a); r2=norm3(b); r=norm3(c)
  th1=acos(dot_product(a,c)/(r1*r))
  th2=acos(dot_product(b,c)/(r2*r))
  tau=acos(dot_product(cross3(a,c),cross3(b,c))/ &
           (norm3(cross3(a,c))*norm3(cross3(b,c))))
  call zpot(r1,r2,r,th1,th2,tau,value)"""
    elif plan.strategy == "multistate-cartesian":
        xx_assignment = (
            f"xx=transpose(coords)*{BOHR_TO_ANGSTROM}"
            if plan.record["id"] == "ch3nh2"
            else "xx=transpose(coords)"
        )
        body = f"""\
  type(dual2) :: xx(3,natoms),uu(5,5),guu(3,natoms,5,5)
  type(dual2) :: vv(5),gvv(3,natoms,5),dvec(3,natoms,5,5),cc(5,5)
  {xx_assignment}
  call zpot(0,xx,uu,guu,vv,gvv,dvec,cc,0)
  value=vv(1)"""
    elif plan.strategy == "nh3-cartesian":
        body = """\
  type(dual2) :: flat(3*natoms),u11,u22,u12,v1,v2
  type(dual2) :: gu11(3*natoms),gu22(3*natoms),gu12(3*natoms)
  type(dual2) :: gv1(3*natoms),gv2(3*natoms)
  flat=reshape(transpose(coords),[3*natoms])
  call zpot(flat,u11,u22,u12,v1,v2,gu11,gu22,gu12,gv1,gv2)
  value=v1"""
    else:
        raise AssertionError(plan.strategy)

    rendered_body = body
    return f"""\
subroutine potlibc_generated_energy(coords,value,natoms)
  use dual2_mod
{module_use}  implicit none
  integer, intent(in) :: natoms
  type(dual2), intent(in) :: coords(natoms,3)
  type(dual2), intent(out) :: value
  interface
    function distance3(left,right) result(distance)
      use dual2_mod
      type(dual2), intent(in) :: left(3),right(3)
      type(dual2) :: distance
    end function distance3
    function norm3(vector) result(length)
      use dual2_mod
      type(dual2), intent(in) :: vector(3)
      type(dual2) :: length
    end function norm3
    function cross3(left,right) result(product)
      use dual2_mod
      type(dual2), intent(in) :: left(3),right(3)
      type(dual2) :: product(3)
    end function cross3
    subroutine pair_distances(coords,natoms,r)
      use dual2_mod
      integer, intent(in) :: natoms
      type(dual2), intent(in) :: coords(natoms,3)
      type(dual2), intent(out) :: r(natoms*(natoms-1)/2)
    end subroutine pair_distances
  end interface
{rendered_body}
end subroutine potlibc_generated_energy

function distance3(left,right) result(distance)
  use dual2_mod
  implicit none
  type(dual2), intent(in) :: left(3),right(3)
  type(dual2) :: distance
  distance=sqrt((left(1)-right(1))**2+(left(2)-right(2))**2+ &
                (left(3)-right(3))**2)
end function distance3

function norm3(vector) result(length)
  use dual2_mod
  implicit none
  type(dual2), intent(in) :: vector(3)
  type(dual2) :: length
  length=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
end function norm3

function cross3(left,right) result(product)
  use dual2_mod
  implicit none
  type(dual2), intent(in) :: left(3),right(3)
  type(dual2) :: product(3)
  product(1)=left(2)*right(3)-left(3)*right(2)
  product(2)=left(3)*right(1)-left(1)*right(3)
  product(3)=left(1)*right(2)-left(2)*right(1)
end function cross3

subroutine pair_distances(coords,natoms,r)
  use dual2_mod
  implicit none
  integer, intent(in) :: natoms
  type(dual2), intent(in) :: coords(natoms,3)
  type(dual2), intent(out) :: r(natoms*(natoms-1)/2)
  integer :: left,right,k
  k=0
  do left=1,natoms-1
    do right=left+1,natoms
      k=k+1
      r(k)=sqrt((coords(left,1)-coords(right,1))**2+ &
                (coords(left,2)-coords(right,2))**2+ &
                (coords(left,3)-coords(right,3))**2)
    end do
  end do
end subroutine pair_distances
"""


def _cbe_adapter() -> str:
    return """\
subroutine potlib_plugin_initialize() bind(c)
  implicit none
  call prepot()
  call zprepot()
end subroutine potlib_plugin_initialize

integer(c_int) function potlib_plugin_atom_count() bind(c)
  use iso_c_binding, only: c_int
  implicit none
  potlib_plugin_atom_count=6
end function potlib_plugin_atom_count

integer(c_int) function potlib_plugin_metadata(natoms,derivatives,atomic_numbers,reference) bind(c)
  use iso_c_binding, only: c_double,c_int
  implicit none
  integer(c_int), value :: natoms
  integer(c_int), intent(out) :: derivatives,atomic_numbers(*)
  real(c_double), intent(out) :: reference(*)
  integer, parameter :: maxatom=25
  character(75) :: ref(5)
  double precision :: cartnu(maxatom,3)
  integer :: indexes(maxatom),irctnt,native_atoms,icartr,mder,msurf
  integer :: atom,axis,k
  common /infocm/ cartnu,indexes,irctnt,native_atoms,icartr,mder,msurf,ref
  if (natoms.ne.native_atoms) then
    potlib_plugin_metadata=1
    return
  end if
  derivatives=2
  do atom=1,natoms
    atomic_numbers(atom)=indexes(atom)
    do axis=1,3
      k=3*(atom-1)+axis
      reference(k)=cartnu(atom,axis)
    end do
  end do
  potlib_plugin_metadata=0
end function potlib_plugin_metadata

subroutine potlib_plugin_evaluate(n,x,energy,force) bind(c)
  use iso_c_binding, only: c_double,c_int
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*)
  real(c_double), intent(out) :: energy(*),force(*)
  integer, parameter :: maxatom=25,isurf=5,jsurf=15,natoms=6
  double precision :: pengygs,pengyes(isurf),pengyij(jsurf)
  double precision :: dgscart(maxatom,3),descart(maxatom,3,isurf)
  double precision :: dijcart(maxatom,3,jsurf),cart(maxatom,3),anuzero
  integer :: nullbl(maxatom),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: sample,atom,axis,k
  common /usrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /usricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+3*(atom-1)+axis
        cart(atom,axis)=x(k)
      end do
    end do
    call pot()
    energy(sample)=pengygs
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+3*(atom-1)+axis
        force(k)=-dgscart(atom,axis)
      end do
    end do
  end do
end subroutine potlib_plugin_evaluate

subroutine potlib_plugin_energy_directional(n,x,direction,derivative) bind(c)
  use dual2_mod
  use iso_c_binding, only: c_double,c_int
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*),direction(*)
  real(c_double), intent(out) :: derivative(*)
  integer, parameter :: maxatom=25,isurf=5,jsurf=15,natoms=6
  type(dual2) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual2) :: dgscart(maxatom,3),descart(maxatom,3,isurf)
  type(dual2) :: dijcart(maxatom,3,jsurf),cart(maxatom,3),anuzero
  integer :: nullbl(maxatom),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: sample,atom,axis,k
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+3*(atom-1)+axis
        cart(atom,axis)=dual2(x(k),0.0d0,direction(k),0.0d0)
      end do
    end do
    call zpot()
    derivative(sample)=pengygs%b
  end do
end subroutine potlib_plugin_energy_directional

subroutine potlib_plugin_hvp(n,x,direction,product) bind(c)
  use dual2_mod
  use iso_c_binding, only: c_double,c_int
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*),direction(*)
  real(c_double), intent(out) :: product(*)
  integer, parameter :: maxatom=25,isurf=5,jsurf=15,natoms=6
  type(dual2) :: pengygs,pengyes(isurf),pengyij(jsurf)
  type(dual2) :: dgscart(maxatom,3),descart(maxatom,3,isurf)
  type(dual2) :: dijcart(maxatom,3,jsurf),cart(maxatom,3),anuzero
  integer :: nullbl(maxatom),nflag(20),nasurf(isurf+1,isurf+1),nder
  integer :: sample,atom,axis,k
  common /zusrocm/ pengygs,pengyes,pengyij,dgscart,descart,dijcart
  common /zusricm/ cart,anuzero,nullbl,nflag,nasurf,nder
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+3*(atom-1)+axis
        cart(atom,axis)=dual2(x(k),0.0d0,direction(k),0.0d0)
      end do
    end do
    call zpot()
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+3*(atom-1)+axis
        product(k)=dgscart(atom,axis)%b
      end do
    end do
  end do
end subroutine potlib_plugin_hvp
"""


def _adapter(plan: SurfacePlan) -> str:
    if plan.record["path"] == "src/CH5/ch52008.f":
        return _cbe_adapter()
    fixed_atoms = 0 if plan.variable_atoms else plan.atoms
    init_name = _initialization_name(plan)
    if plan.strategy in {"triatomic-function", "distance-function"}:
        init_name = "potlibc_standalone_init"
        init_use = (
            "  use potlibc_standalone_pes, only: "
            "zpotlibc_standalone_init\n"
        )
    else:
        init_use = (
            f"  use {plan.module_name}, only: z{init_name}\n"
            if plan.module_name and init_name
            else ""
        )
    init_call = f"  call z{init_name}()\n" if init_name else ""
    metadata_init = ""
    if plan.record["id"] == "n4pes-l-imls-g2-v1":
        init_use = (
            "  use potlibc_n4_data_state, only: zpotlibc_n4_prepare\n"
        )
        init_call = "  call zpotlibc_n4_prepare()\n  call zprepot()\n"
    if plan.record["id"] in {"nich4", "niheam5", "niheam6", "niheam7"}:
        archived_units = (
            (4, 7, 9, 10, 93, 94, 95, 96)
            if plan.record["id"] == "nich4"
            else (4, 9)
        )
        init_call = "".join(
            f"  open(unit={unit},status='scratch',action='readwrite')\n"
            f"  call zpotlibc_archived_unit{unit}({unit})\n"
            f"  rewind({unit})\n"
            for unit in archived_units
        ) + init_call
        metadata_init = "  call zsetup(3*natoms)\n"
    atom_check = (
        "  if (natoms.lt.1) then"
        if plan.variable_atoms
        else f"  if (natoms.ne.{fixed_atoms}) then"
    )
    metadata = plan.record.get("atomic_metadata", {"kind": "unknown"})
    fixed_numbers = (
        metadata.get("numbers")
        if metadata.get("kind") == "fixed" and isinstance(metadata.get("numbers"), list)
        else None
    )
    variable_number = (
        int(metadata["number"])
        if metadata.get("kind") == "homonuclear-variable" and "number" in metadata
        else None
    )
    if fixed_numbers is not None and len(fixed_numbers) == fixed_atoms and fixed_atoms > 0:
        atom_metadata = "\n".join(
            f"  atomic_numbers({index})={number}"
            for index, number in enumerate(fixed_numbers, 1)
        )
    elif variable_number is not None:
        atom_metadata = f"""\
  do index=1,natoms
    atomic_numbers(index)={variable_number}
  end do"""
    elif plan.record["id"] in {"vbo-beta1", "vbo2-beta1"}:
        atom_metadata = """\
  do index=1,natoms
    atomic_numbers(index)=13
  end do
"""
    elif plan.record["id"] == "sioh":
        atom_metadata = """\
  do index=1,natoms
    select case (mod(index-1,3))
    case (0)
      atomic_numbers(index)=14
    case (1)
      atomic_numbers(index)=8
    case (2)
      atomic_numbers(index)=1
    end select
  end do"""
    else:
        atom_metadata = """\
  do index=1,natoms
    atomic_numbers(index)=0
  end do"""
    return f"""\
module potlibc_generated_state
  use iso_c_binding, only: c_int
  implicit none
  integer(c_int), save :: selected_atoms={fixed_atoms}
end module potlibc_generated_state

subroutine potlib_plugin_initialize() bind(c)
{init_use}  implicit none
{init_call}end subroutine potlib_plugin_initialize

integer(c_int) function potlib_plugin_atom_count() bind(c)
  use iso_c_binding, only: c_int
  implicit none
  potlib_plugin_atom_count={fixed_atoms}
end function potlib_plugin_atom_count

integer(c_int) function potlib_plugin_metadata(natoms,derivatives,atomic_numbers,reference) bind(c)
  use iso_c_binding, only: c_double,c_int
  use potlibc_generated_state, only: selected_atoms
  implicit none
  integer(c_int), value :: natoms
  integer(c_int), intent(out) :: derivatives,atomic_numbers(*)
  real(c_double), intent(out) :: reference(*)
  integer :: index
  interface
    subroutine potlibc_state_checkpoint() bind(c)
    end subroutine potlibc_state_checkpoint
  end interface
{atom_check}
    potlib_plugin_metadata=1
    return
  end if
  selected_atoms=natoms
{metadata_init}  derivatives=2
{atom_metadata}
  do index=1,3*natoms
    reference(index)=0.0_c_double
  end do
  call potlibc_state_checkpoint()
  potlib_plugin_metadata=0
end function potlib_plugin_metadata

subroutine potlib_plugin_evaluate(n,x,energy,force) bind(c)
  use dual2_mod
  use iso_c_binding, only: c_double,c_int
  use potlibc_generated_state, only: selected_atoms
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*)
  real(c_double), intent(out) :: energy(*),force(*)
  type(dual2), allocatable :: coords(:,:)
  type(dual2) :: value
  integer :: sample,target,atom,axis,k
  integer(c_int) :: natoms
  natoms=selected_atoms
  allocate(coords(natoms,3))
  do sample=1,n
    do target=1,3*natoms
      do atom=1,natoms
        do axis=1,3
          k=(sample-1)*3*natoms+3*(atom-1)+axis
          coords(atom,axis)%re=x(k)
          coords(atom,axis)%a=0.0d0
          coords(atom,axis)%b=merge(1.0d0,0.0d0,3*(atom-1)+axis==target)
          coords(atom,axis)%ab=0.0d0
        end do
      end do
      call potlibc_generated_restore()
      call potlibc_generated_energy(coords,value,natoms)
      energy(sample)=value%re
      k=(sample-1)*3*natoms+target
      force(k)=-value%b
    end do
  end do
end subroutine potlib_plugin_evaluate

subroutine potlib_plugin_energy_directional(n,x,direction,derivative) bind(c)
  use dual2_mod
  use iso_c_binding, only: c_double,c_int
  use potlibc_generated_state, only: selected_atoms
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*),direction(*)
  real(c_double), intent(out) :: derivative(*)
  type(dual2), allocatable :: coords(:,:)
  type(dual2) :: value
  integer :: sample,atom,axis,k
  integer(c_int) :: natoms
  natoms=selected_atoms
  allocate(coords(natoms,3))
  do sample=1,n
    do atom=1,natoms
      do axis=1,3
        k=(sample-1)*3*natoms+3*(atom-1)+axis
        coords(atom,axis)%re=x(k)
        coords(atom,axis)%a=0.0d0
        coords(atom,axis)%b=direction(k)
        coords(atom,axis)%ab=0.0d0
      end do
    end do
    call potlibc_generated_restore()
    call potlibc_generated_energy(coords,value,natoms)
    derivative(sample)=value%b
  end do
end subroutine potlib_plugin_energy_directional

subroutine potlib_plugin_hvp(n,x,direction,product) bind(c)
  use dual2_mod
  use iso_c_binding, only: c_double,c_int
  use potlibc_generated_state, only: selected_atoms
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: x(*),direction(*)
  real(c_double), intent(out) :: product(*)
  type(dual2), allocatable :: coords(:,:)
  type(dual2) :: value
  integer :: sample,target,atom,axis,k
  integer(c_int) :: natoms
  natoms=selected_atoms
  allocate(coords(natoms,3))
  do sample=1,n
    do target=1,3*natoms
      do atom=1,natoms
        do axis=1,3
          k=(sample-1)*3*natoms+3*(atom-1)+axis
          coords(atom,axis)%re=x(k)
          coords(atom,axis)%a=merge(1.0d0,0.0d0,3*(atom-1)+axis==target)
          coords(atom,axis)%b=direction(k)
          coords(atom,axis)%ab=0.0d0
        end do
      end do
      call potlibc_generated_restore()
      call potlibc_generated_energy(coords,value,natoms)
      k=(sample-1)*3*natoms+target
      product(k)=value%ab
    end do
  end do
end subroutine potlib_plugin_hvp

subroutine potlibc_generated_restore()
  implicit none
  interface
    subroutine potlibc_state_restore() bind(c)
    end subroutine potlibc_state_restore
  end interface
  call potlibc_state_restore()
end subroutine potlibc_generated_restore

{_energy_bridge(plan)}"""


def generate(
    vendor: Path,
    surface_id: str,
    output_source: Path,
    output_utility: Path,
    output_adapter: Path,
    output_data_c: Path | None = None,
    output_data_bin: Path | None = None,
) -> SurfacePlan:
    catalog = inventory(vendor)
    matches = [item for item in catalog["surfaces"] if item["id"] == surface_id]
    if not matches:
        raise ValueError(f"unknown surface id {surface_id!r}")
    plan = _plan(matches[0], vendor)
    transformed_source, transformed_utility = _transform(plan, vendor)
    for output in (output_source, output_utility, output_adapter):
        output.parent.mkdir(parents=True, exist_ok=True)
    output_source.write_text(transformed_source, encoding="utf-8")
    output_utility.write_text(transformed_utility, encoding="utf-8")
    output_adapter.write_text(_adapter(plan), encoding="utf-8")
    if output_data_c is not None and output_data_bin is not None:
        output_data_c.parent.mkdir(parents=True, exist_ok=True)
        output_data_bin.parent.mkdir(parents=True, exist_ok=True)
        if surface_id == "n4pes-l-imls-g2-v1":
            package_path = plan.source.with_suffix(".tar")
            member = f"{package_path.stem}/pes.dat"
            with tarfile.open(package_path) as package:
                extracted = package.extractfile(member)
                if extracted is None:
                    raise ValueError(f"{surface_id}: archive has no {member}")
                output_data_bin.write_bytes(extracted.read())
            binary_path = str(output_data_bin.resolve()).replace("\\", "\\\\")
            output_data_c.write_text(
                f"""\
#define _GNU_SOURCE
#include <errno.h>
#include <linux/memfd.h>
#include <stddef.h>
#include <sys/mman.h>
#include <unistd.h>

__asm__(".pushsection .rodata\\n"
        ".global potlibc_n4_data_start\\n"
        "potlibc_n4_data_start:\\n"
        ".incbin \\"{binary_path}\\"\\n"
        ".global potlibc_n4_data_end\\n"
        "potlibc_n4_data_end:\\n"
        ".popsection\\n");
extern const unsigned char potlibc_n4_data_start[];
extern const unsigned char potlibc_n4_data_end[];
static int data_fd = -1;

int potlibc_n4_data_fd(void)
{{
    if (data_fd >= 0) {{
        return data_fd;
    }}
    data_fd = memfd_create("potlibc-n4-data", MFD_CLOEXEC);
    if (data_fd < 0) {{
        return -1;
    }}
    const unsigned char *position = potlibc_n4_data_start;
    size_t remaining =
        (size_t)(potlibc_n4_data_end - potlibc_n4_data_start);
    while (remaining != 0) {{
        ssize_t written = write(data_fd, position, remaining);
        if (written < 0 && errno == EINTR) {{
            continue;
        }}
        if (written <= 0) {{
            close(data_fd);
            data_fd = -1;
            return -1;
        }}
        position += written;
        remaining -= (size_t)written;
    }}
    if (lseek(data_fd, 0, SEEK_SET) < 0) {{
        close(data_fd);
        data_fd = -1;
    }}
    return data_fd;
}}

__attribute__((destructor)) static void close_data(void)
{{
    if (data_fd >= 0) close(data_fd);
}}
""",
                encoding="utf-8",
            )
        else:
            output_data_c.write_text(
                "/* No embedded binary data for this surface. */\n",
                encoding="utf-8",
            )
            output_data_bin.write_bytes(b"")
    return plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a mixed-forward-dual POTLIB plugin at build time."
    )
    parser.add_argument("vendor_root", type=Path)
    parser.add_argument("surface_id")
    parser.add_argument("output_source", type=Path)
    parser.add_argument("output_utility", type=Path)
    parser.add_argument("output_adapter", type=Path)
    parser.add_argument("output_data_c", type=Path, nargs="?")
    parser.add_argument("output_data_bin", type=Path, nargs="?")
    args = parser.parse_args(argv)
    try:
        plan = generate(
            args.vendor_root,
            args.surface_id,
            args.output_source,
            args.output_utility,
            args.output_adapter,
            args.output_data_c,
            args.output_data_bin,
        )
    except (OSError, RuntimeError, ValueError, StopIteration) as error:
        parser.error(str(error))
    print(
        f"{plan.record['id']}: {plan.record['interface_family']} "
        f"via {plan.strategy}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
