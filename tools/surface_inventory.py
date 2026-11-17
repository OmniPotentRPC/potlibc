#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Inventory every Fortran source in the vendored PotLib 2023 archive."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True

FORTRAN_SUFFIXES = {".f", ".F", ".f90", ".F90"}
EXPECTED_FILE_COUNT = 403
HEADER = re.compile(r"(?im)^[c*!][ \t]*([^:\n]+?)[ \t]*:[ \t]*(.*?)[ \t]*$")
DOI = re.compile(
    r"(?i)\b10\.\d{4,9}/[-._;()/:A-Z0-9]+"
)
ROUTINE = re.compile(
    r"(?i)^\s*(?:[a-z][\w*() ]+\s+)?"
    r"(subroutine|function|program|module)\s+([a-z_]\w*)\s*(?:\((.*)\))?\s*$"
)
ENTRY = re.compile(r"(?i)^\s*entry\s+([a-z_]\w*)\s*(?:\((.*)\))?\s*$")
USE = re.compile(
    r"(?im)^\s*use(?:\s*,[^:]+::)?\s+([a-z_]\w*)", re.IGNORECASE
)
INCLUDE = re.compile(r"(?im)^\s*include\s*[\"']([^\"']+)[\"']")
CALL = re.compile(r"(?i)\bcall\s+([a-z_]\w*)")
COMMON = re.compile(r"(?i)\bcommon\s*/\s*([a-z_]\w*)\s*/")
FILE_ARGUMENT = re.compile(r"(?i)\bfile\s*=\s*(['\"])(.*?)\1")
UNIT_ARGUMENT = re.compile(r"(?i)(?:\bunit\s*=\s*)?(\d+)")
READ = re.compile(r"(?i)\bread\s*\(\s*(?:unit\s*=\s*)?(\d+)")

KNOWN_HELPERS = {
    "C2H7/interface_c2h7pot5_mod.f90": "adapter_module",
    "C3H7NO/Full-PES/basis.f90": "helper_module",
    "C3H7NO/Full-PES/gradient.f90": "helper_module",
    "CH4/getpot_ch4t8.f90": "adapter",
    "done/FH2O/jcb3_xyz.f": "coordinate_helper",
}
KNOWN_CUSTOM_SURFACES = {
    "F2H2/bjkkl.f",
    "done/FH2O/soc-fh2o-nn.f",
}
ARCHIVED_SURFACES = {
    "NinR/nich4.f": ("NinR/nich4.tar", "METHANE/methane.f"),
    "NinR/niheam5.f": ("NinR/niheam5.tar", "EAM5/hni5.f"),
    "NinR/niheam6.f": ("NinR/niheam6.tar", "EAM6/hni6.f"),
    "NinR/niheam7.f": ("NinR/niheam7.tar", "EAM7/hni7.f"),
}

ATOMIC_NUMBERS = {
    "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "NE": 10, "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15,
    "S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21, "TI": 22,
    "V": 23, "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29,
    "ZN": 30, "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36,
    "RB": 37, "SR": 38, "Y": 39, "ZR": 40, "NB": 41, "MO": 42, "TC": 43,
    "RU": 44, "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, "SN": 50,
    "SB": 51, "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "LA": 57,
    "CE": 58, "PR": 59, "ND": 60, "PM": 61, "SM": 62, "EU": 63, "GD": 64,
    "TB": 65, "DY": 66, "HO": 67, "ER": 68, "TM": 69, "YB": 70, "LU": 71,
    "HF": 72, "TA": 73, "W": 74, "RE": 75, "OS": 76, "IR": 77, "PT": 78,
    "AU": 79, "HG": 80, "TL": 81, "PB": 82, "BI": 83, "PO": 84, "AT": 85,
    "RN": 86,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_source(path: Path) -> str:
    return path.read_bytes().decode("latin-1")


def header_fields(text: str) -> dict[str, str]:
    """Return normalized PotLib header fields from the beginning of a file."""
    fields: dict[str, str] = {}
    for key, value in HEADER.findall(text[:16000]):
        key = " ".join(key.lower().split())
        value = " ".join(value.split())
        if value:
            fields[key] = value
    return fields


def header_block(text: str, names: tuple[str, ...]) -> str:
    """Return a normalized, possibly multiline field from a PotLib header."""
    wanted = {name.lower() for name in names}
    values: list[str] = []
    active = False
    for raw in text[:16000].splitlines():
        if not raw or raw[0] not in "cC*!":
            if active and raw.strip():
                break
            continue
        body = raw[1:]
        line = body.strip()
        match = re.match(r"([^:]+?)\s*:+\s*(.*)$", line)
        if match and len(body) - len(body.lstrip()) <= 4:
            key = " ".join(match.group(1).lower().split())
            active = key in wanted
            if active and match.group(2).strip():
                values.append(match.group(2).strip())
            continue
        if active:
            if not line:
                active = False
            else:
                values.append(line)
    return " ".join(" ".join(values).split())


def surface_documentation(text: str, fields: dict[str, str]) -> dict[str, Any]:
    """Extract the deposited descriptive metadata without interpreting it."""
    reference_parts = [
        header_block(text, ("reference", "references", "ref")),
        header_block(text, ("cross reference",)),
    ]
    reference = "; ".join(part for part in reference_parts if part)
    dois = sorted(
        {
            match.group(0).rstrip(".,;)")
            for match in DOI.finditer(text[:16000])
        }
    )
    return {
        "formula": fields.get("system", ""),
        "coordinates": header_block(
            text,
            (
                "coordinates",
                "coordinates - distances",
                "coord",
                "input coordinate matrix",
            ),
        ),
        "reference": reference,
        "dois": dois,
    }


def active_lines(text: str, fixed_form: bool) -> list[str]:
    """Remove fixed/free-form comments and join Fortran continuations."""
    physical: list[str] = []
    for raw in text.splitlines():
        if not raw.strip():
            continue
        if raw[0] in "cC*":
            continue
        stripped = raw.lstrip()
        if stripped.startswith("!"):
            continue
        physical.append(raw.split("!", 1)[0].rstrip())

    fixed: list[str] = []
    for raw in physical:
        fixed_continuation = fixed_form and (
            len(raw) > 5 and raw[5] not in {" ", "0"} and not raw[:5].strip()
        )
        if fixed_continuation and fixed:
            fixed[-1] += " " + raw[6:].strip()
        else:
            fixed.append(raw.strip())

    logical: list[str] = []
    pending = ""
    for raw in fixed:
        free_continuation = raw.lstrip().startswith("&")
        body = raw.lstrip()[1:].lstrip() if free_continuation else raw.strip()
        if pending:
            pending += " " + body
        else:
            pending = body
        if pending.endswith("&"):
            pending = pending[:-1].rstrip()
        else:
            logical.append(" ".join(pending.split()))
            pending = ""
    if pending:
        logical.append(" ".join(pending.split()))
    return logical


def declarations(lines: list[str]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for line in lines:
        match = ENTRY.match(line)
        if match:
            result.append(
                {
                    "kind": "entry",
                    "name": match.group(1),
                    "arguments": split_arguments(match.group(2)),
                    "signature": line,
                }
            )
            continue
        match = ROUTINE.match(line)
        if match:
            result.append(
                {
                    "kind": match.group(1).lower(),
                    "name": match.group(2),
                    "arguments": split_arguments(match.group(3)),
                    "signature": line,
                }
            )
    return result


def split_arguments(arguments: str | None) -> list[str]:
    if not arguments:
        return []
    return [part.strip() for part in arguments.split(",") if part.strip()]


def stable_id(relative: str) -> str:
    stem = Path(relative).stem.lower()
    slug = re.sub(r"[^a-z0-9]+", "-", stem).strip("-")
    if relative == "N4/n4pes-gpip-v2/n4pes-gpip-v2.f":
        return "n4pes-gpip-v2-package"
    return slug


def interface_documents(vendor: Path) -> dict[str, dict[str, str]]:
    documents: dict[str, dict[str, str]] = {}
    for path in sorted((vendor / "interface").glob("*.html")):
        text = path.read_text(encoding="latin-1", errors="replace")
        heading = re.search(r"(?is)<h2[^>]*>\s*(?:Interface\s+)?(.*?)\s*</h2>", text)
        pre = re.search(r"(?is)<pre[^>]*>(.*?)</pre>", text)
        if not heading:
            continue
        name = re.sub(r"<[^>]+>", "", heading.group(1)).strip()
        family = next(
            (
                candidate
                for candidate in (
                    "POTLIB-2001",
                    "HO-MM-0",
                    "HO-MM-1",
                    "3-1V",
                    "3-2V",
                    "3-XS",
                    "4-1S",
                    "4-XS",
                    "5-1S",
                    "6-1S",
                    "7-1S",
                )
                if candidate.lower() in name.lower()
            ),
            name,
        )
        sample = ""
        if pre:
            sample = " ".join(re.sub(r"<[^>]+>", "", pre.group(1)).split())
        documents[family.upper()] = {
            "path": path.relative_to(vendor).as_posix(),
            "sample_signature": sample,
        }
    return documents


def choose_entry(
    items: list[dict[str, Any]], preferred_family: str = ""
) -> dict[str, Any] | None:
    priorities = (
        ("surf", "pot")
        if preferred_family.upper().replace(" ", "-") == "SECTION-2"
        else ("pot", "surf")
    ) + (
        "getpot",
        "pesh3cl",
        "ho3_pes_interface",
        "n4pipnn",
        "fh2onn",
        "h3o2pipnn",
        "pes",
        "f",
    )
    lowered = [
        (item["name"].lower(), item)
        for item in items
        if item["kind"] in {"subroutine", "function", "entry"}
    ]
    for name in priorities:
        for candidate, item in lowered:
            if candidate == name:
                return item
    for candidate, item in lowered:
        if any(token in candidate for token in ("pot", "pes", "energy_gradient")):
            return item
    return lowered[0][1] if lowered else None


def interface_family(
    fields: dict[str, str],
    entry: dict[str, Any] | None,
    commons: set[str],
    is_module_package: bool,
) -> tuple[str, str]:
    declared = fields.get("interface", "").strip()
    if declared:
        normalized = declared.upper().replace("SECTION 2", "SECTION-2")
        normalized = {
            "POTLIB200": "POTLIB-2001",
            "POTLIB2001": "POTLIB-2001",
        }.get(normalized, normalized)
        return normalized, "header"
    if entry:
        name = entry["name"].lower()
        args = [arg.lower() for arg in entry["arguments"]]
        joined = ",".join(args)
        if name == "surf":
            return "SECTION-2", "entry-name"
        if {"natom", "maxatom"}.issubset(args) and all(
            coordinate in args for coordinate in ("x", "y", "z")
        ):
            derivative = any(arg in joined for arg in ("der", "grad", "dvd"))
            return ("HO-MM-1" if derivative else "HO-MM-0"), "signature"
        if name == "pot" and {"r", "e", "nsurf"} <= set(args) and any(
            count in args for count in ("nt", "nvals")
        ):
            derivative = any(arg in joined for arg in ("de", "dedr", "grad"))
            return ("3-2V" if derivative else "3-1V"), "signature"
        if name == "pot" and args[:2] in (["rvp", "evp"], ["rr", "ee"]):
            return "3-1V", "signature"
    lowered_commons = {name.lower() for name in commons}
    if {"infocm", "usricm"} <= lowered_commons:
        return "POTLIB-2001", "common-blocks"
    if entry and entry["name"].lower() == "pot" and not entry["arguments"]:
        if "pt1cm" in lowered_commons:
            return "3-XS", "common-blocks"
        if "potcm" in lowered_commons:
            return "3-1S", "common-blocks"
    if is_module_package:
        return "CUSTOM-PACKAGE", "module/package-entry"
    return "SECTION-2", "nonstandard-entry"


def atom_count(
    fields: dict[str, str],
    entry: dict[str, Any] | None,
    family: str,
    relative: str,
) -> dict[str, Any]:
    value = fields.get("number of bodies", "")
    match = re.search(r"\d+", value)
    if match:
        return {"kind": "fixed", "value": int(match.group())}
    if "variable" in value.lower():
        return {"kind": "variable", "symbol": variable_atom_symbol(entry)}
    family_atoms = re.match(r"(\d+)-", family)
    if family_atoms:
        return {"kind": "fixed", "value": int(family_atoms.group(1))}
    if entry:
        for argument in entry["arguments"]:
            if argument.lower() in {"natom", "natoms"}:
                return {"kind": "variable", "symbol": argument.upper()}
        signature = entry["signature"]
        dimensions = re.findall(r"(?i)\(\s*3\s*,\s*(\d+)\s*\)", signature)
        if dimensions:
            return {"kind": "fixed", "value": int(dimensions[0])}
    for candidate in (fields.get("system", ""), package_key(relative).split("/")[-1]):
        inferred = formula_atom_count(candidate)
        if inferred is not None:
            return {"kind": "fixed", "value": inferred, "inferred_from": candidate}
    return {"kind": "unknown", "value": None}


def variable_atom_symbol(entry: dict[str, Any] | None) -> str:
    if entry:
        for argument in entry["arguments"]:
            if argument.lower() in {"natom", "natoms"}:
                return argument.upper()
    return "NATOM"


def call_argument_lists(text: str, name: str) -> list[str]:
    results: list[str] = []
    for match in re.finditer(rf"(?i)\b{re.escape(name)}\s*\(", text):
        start = match.end()
        depth = 1
        quote = ""
        for index in range(start, len(text)):
            character = text[index]
            if quote:
                if character == quote:
                    quote = ""
                continue
            if character in {"'", '"'}:
                quote = character
            elif character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    results.append(text[start:index])
                    break
    return results


def keyword_expression(arguments: str, keyword: str) -> str | None:
    match = re.search(rf"(?i)\b{re.escape(keyword)}\s*=\s*", arguments)
    if not match:
        return None
    start = match.end()
    depth = 0
    quote = ""
    for index in range(start, len(arguments)):
        character = arguments[index]
        if quote:
            if character == quote:
                quote = ""
            continue
        if character in {"'", '"'}:
            quote = character
        elif character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        elif character == "," and depth == 0:
            return arguments[start:index].strip()
    return arguments[start:].strip()


def formula_atom_count(system: str) -> int | None:
    compact = re.sub(r"[^A-Za-z0-9+]", "", system)
    if not compact or len(compact) > 24:
        return None
    if compact.isalpha() and compact.isupper():
        return len(compact)
    parts = compact.split("+")
    total = 0
    for part in parts:
        position = 0
        subtotal = 0
        for match in re.finditer(r"([A-Z][a-z]?)(\d*)", part):
            if match.start() != position:
                return None
            subtotal += int(match.group(2) or "1")
            position = match.end()
        if position != len(part) or subtotal == 0:
            return None
        total += subtotal
    return total or None


def formula_symbols(system: str) -> list[str] | None:
    compact = re.sub(r"[^A-Za-z0-9+]", "", system)
    if not compact or len(compact) > 32:
        return None
    symbols: list[str] = []
    for part in compact.split("+"):
        position = 0
        for match in re.finditer(r"([A-Z][a-z]?)(\d*)", part):
            if match.start() != position:
                return None
            symbol = match.group(1).upper()
            count = int(match.group(2) or "1")
            symbols.extend([symbol] * count)
            position = match.end()
        if position != len(part):
            return None
    return symbols or None


def atomic_metadata(
    system: str, count: dict[str, Any]
) -> dict[str, Any]:
    symbols = formula_symbols(system)
    if count["kind"] == "fixed":
        expected = int(count["value"])
        if symbols is not None and len(symbols) == expected:
            if all(symbol in ATOMIC_NUMBERS for symbol in symbols):
                return {
                    "kind": "fixed",
                    "source": "system_formula",
                    "symbols": symbols,
                    "numbers": [ATOMIC_NUMBERS[symbol] for symbol in symbols],
                }
            return {
                "kind": "unknown",
                "reason": "formula contains symbols outside periodic table map",
            }
        return {
            "kind": "unknown",
            "reason": "system formula does not resolve a unique fixed atom list",
        }
    if count["kind"] == "variable":
        if symbols is not None and len(set(symbols)) == 1 and symbols[0] in ATOMIC_NUMBERS:
            return {
                "kind": "homonuclear-variable",
                "source": "system_formula",
                "symbol": symbols[0],
                "number": ATOMIC_NUMBERS[symbols[0]],
            }
        return {
            "kind": "unknown",
            "reason": "variable-size system has no unambiguous elemental species",
        }
    return {"kind": "unknown", "reason": "atom count is unknown"}


def external_data_references(
    text: str, fields: dict[str, str], source: Path, vendor: Path
) -> tuple[list[dict[str, Any]], list[int]]:
    active = "\n".join(active_lines(text, source.suffix.lower() != ".f90"))
    read_units = {int(unit) for unit in READ.findall(active)}
    references: dict[str, dict[str, Any]] = {}
    declared = fields.get("data file", "")
    if declared:
        for token in re.findall(r"[\w.+-]+\.(?:dat|txt|par|tar|in|pes)\b", declared, re.I):
            references[token] = {"name": token, "evidence": ["header"]}
    for contents in call_argument_lists(active, "open"):
        expression = keyword_expression(contents, "file")
        if not expression:
            continue
        file_match = FILE_ARGUMENT.fullmatch(f"file={expression}")
        literal = file_match is not None
        name = file_match.group(2) if file_match else expression
        record = references.setdefault(name, {"name": name, "evidence": []})
        record["evidence"].append("OPEN")
        record["literal"] = literal
        unit_match = UNIT_ARGUMENT.search(contents)
        if unit_match:
            unit = int(unit_match.group(1))
            record["unit"] = unit
            record["read"] = unit in read_units
    for record in references.values():
        matches = []
        if record.get("literal", True):
            matches = sorted(
                path.relative_to(vendor).as_posix()
                for path in source.parent.rglob(Path(record["name"]).name)
                if path.is_file()
            )
        record["resolved_paths"] = matches
        record["evidence"] = sorted(set(record["evidence"]))
    return sorted(references.values(), key=lambda item: item["name"].lower()), sorted(read_units)


def package_key(relative: str) -> str:
    parts = relative.split("/")
    if parts[0] == "done" and len(parts) > 1:
        return "/".join(parts[:2])
    return parts[0]


def is_test_driver(relative: str, _decls: list[dict[str, Any]]) -> bool:
    basename = Path(relative).name.lower()
    return basename in {"test.f", "test.f90", "main.f", "main.f90"}


def classify_source(
    relative: str,
    text: str,
    fields: dict[str, str],
    decls: list[dict[str, Any]],
) -> tuple[bool, str, list[str]]:
    evidence: list[str] = []
    if "__MACOSX/" in relative or Path(relative).name.startswith("._"):
        return False, "macos_metadata", evidence
    if relative in KNOWN_HELPERS:
        return False, KNOWN_HELPERS[relative], evidence
    if relative in KNOWN_CUSTOM_SURFACES:
        return True, "package_surface", ["known custom package PES entry"]

    names = {item["name"].lower() for item in decls}
    canonical = "pot" in names
    surf = "surf" in names
    if canonical:
        evidence.append("active POT/ENTRY POT")
    if surf:
        evidence.append("active Section-2 SURF entry")
    if fields.get("interface", "").upper().startswith(("HO-MM", "HE-MM")):
        evidence.append("molecular-mechanics interface header")
    header_surface = bool(fields.get("system") and fields.get("interface"))
    has_implementation = any(
        item["kind"] in {"subroutine", "function"} for item in decls
    )
    package_entry = any(
        any(token in name for token in ("pot", "pes", "energy_gradient"))
        for name in names
    )

    if canonical or surf:
        return True, "surface_entry", evidence
    if header_surface and has_implementation:
        evidence.append("PotLib surface header and active implementation")
        return True, "declared_surface", evidence
    if package_entry and has_implementation:
        evidence.append("custom module/package PES implementation")
        return True, "package_surface", evidence
    if is_test_driver(relative, decls):
        return False, "test_driver", evidence
    if fields.get("interface") and not has_implementation:
        return False, "package_manifest", evidence
    if any(item["kind"] == "module" for item in decls):
        return False, "helper_module", evidence
    if any(item["kind"] == "program" for item in decls):
        if any(item["kind"] == "function" for item in decls):
            evidence.append("standalone package containing PES function")
            return True, "package_surface", evidence
        return False, "test_driver", evidence
    return False, "helper_source", evidence


def build_dependencies(records: dict[str, dict[str, Any]]) -> None:
    modules: dict[str, set[str]] = defaultdict(set)
    symbols: dict[str, set[str]] = defaultdict(set)
    basenames: dict[str, set[str]] = defaultdict(set)
    for relative, record in records.items():
        basenames[Path(relative).name.lower()].add(relative)
        for item in record["_declarations"]:
            if item["kind"] == "module":
                modules[item["name"].lower()].add(relative)
            elif item["kind"] in {"subroutine", "function"}:
                symbols[item["name"].lower()].add(relative)

    for relative, record in records.items():
        text = record["_active"]
        package = package_key(relative)
        dependencies: dict[str, set[str]] = defaultdict(set)
        for name in USE.findall(text):
            for provider in modules.get(name.lower(), set()):
                if provider != relative and package_key(provider) == package:
                    dependencies[provider].add("USE")
        for name in INCLUDE.findall(text):
            candidates = basenames.get(Path(name).name.lower(), set())
            for provider in candidates:
                if provider != relative and package_key(provider) == package:
                    dependencies[provider].add("INCLUDE")
        for name in CALL.findall(text):
            providers = {
                provider
                for provider in symbols.get(name.lower(), set())
                if provider != relative and package_key(provider) == package
            }
            if len(providers) == 1:
                dependencies[next(iter(providers))].add("CALL")
        record["source_dependencies"] = [
            {"path": path, "evidence": sorted(kinds)}
            for path, kinds in sorted(dependencies.items())
        ]


def owner_ids(
    relative: str,
    records: dict[str, dict[str, Any]],
    surface_ids: dict[str, str],
) -> list[str]:
    owners: set[str] = set()
    dependencies = {
        dependency["path"] for dependency in records[relative]["source_dependencies"]
    }
    owners.update(surface_ids[path] for path in dependencies if path in surface_ids)
    for surface_path, identifier in surface_ids.items():
        surface_dependencies = {
            dependency["path"]
            for dependency in records[surface_path]["source_dependencies"]
        }
        if relative in surface_dependencies:
            owners.add(identifier)
    if owners:
        return sorted(owners)

    normalized = relative.replace("__MACOSX/", "").replace("/._", "/")
    parent = Path(normalized).parent
    same_directory = {
        identifier
        for path, identifier in surface_ids.items()
        if Path(path).parent == parent
    }
    if same_directory:
        return sorted(same_directory)
    same_package = {
        identifier
        for path, identifier in surface_ids.items()
        if package_key(path) == package_key(relative)
    }
    return sorted(same_package)


def inventory(vendor: Path) -> dict[str, Any]:
    vendor = vendor.resolve()
    src = vendor / "src"
    if not src.is_dir():
        raise ValueError(f"vendor root has no src directory: {vendor}")
    paths = sorted(
        path
        for path in src.rglob("*")
        if path.is_file() and path.suffix in FORTRAN_SUFFIXES
    )
    documents = interface_documents(vendor)
    records: dict[str, dict[str, Any]] = {}
    ambiguities: list[dict[str, str]] = []

    for path in paths:
        relative = path.relative_to(src).as_posix()
        text = read_source(path)
        fields = header_fields(text)
        implementation_text = text
        archive_record: dict[str, str] | None = None
        if relative in ARCHIVED_SURFACES:
            archive_relative, member = ARCHIVED_SURFACES[relative]
            archive = src / archive_relative
            with tarfile.open(archive) as package:
                extracted = package.extractfile(member)
                if extracted is None:
                    raise ValueError(f"{archive_relative} has no {member}")
                implementation = extracted.read()
            implementation_text = implementation.decode("latin-1")
            archive_record = {
                "path": f"src/{archive_relative}",
                "member": member,
                "archive_sha256": digest(archive),
                "member_sha256": hashlib.sha256(implementation).hexdigest(),
            }
        lines = active_lines(
            implementation_text,
            relative in ARCHIVED_SURFACES or path.suffix.lower() != ".f90",
        )
        decls = declarations(lines)
        if archive_record is not None:
            is_surface = True
            classification = "archived_package_surface"
            evidence = ["implementation source in vendor package archive"]
        else:
            is_surface, classification, evidence = classify_source(
                relative, text, fields, decls
            )
        entry = (
            choose_entry(decls, fields.get("interface", "")) if is_surface else None
        )
        commons = set(COMMON.findall("\n".join(lines)))
        family, family_evidence = interface_family(
            fields,
            entry,
            commons,
            any(item["kind"] == "module" for item in decls),
        )
        external, read_units = external_data_references(
            text, fields, path, vendor
        )
        record: dict[str, Any] = {
            "path": f"src/{relative}",
            "sha256": digest(path),
            "classification": classification,
            "source_dependencies": [],
            "external_data_files": external,
            "_relative": relative,
            "_active": "\n".join(lines),
            "_declarations": decls,
            "_surface": is_surface,
        }
        if is_surface:
            identifier = stable_id(relative)
            record.update(
                {
                    "id": identifier,
                    "system": fields.get("system", ""),
                    "name": fields.get("common name", path.stem),
                    "documentation": surface_documentation(text, fields),
                    "atom_count": atom_count(fields, entry, family, relative),
                    "interface_family": family,
                    "interface_evidence": family_evidence,
                    "interface_document": documents.get(family.upper(), {}).get(
                        "path"
                    ),
                    "entry": entry,
                    "declarations": decls,
                    "common_blocks": sorted(name.lower() for name in commons),
                    "surface_evidence": evidence,
                }
            )
            record["atomic_metadata"] = atomic_metadata(
                record["system"], record["atom_count"]
            )
            if archive_record is not None:
                record["archive_source"] = archive_record
            if record["atom_count"]["kind"] == "unknown":
                ambiguities.append(
                    {
                        "path": record["path"],
                        "code": "unknown_atom_count",
                        "detail": "No body count or variable NATOM signature was found.",
                    }
                )
            if entry is None:
                ambiguities.append(
                    {
                        "path": record["path"],
                        "code": "unknown_entry",
                        "detail": "Surface evidence exists but no callable entry was parsed.",
                    }
                )
        else:
            record.update(
                {
                    "category": classification,
                    "reason": {
                        "macos_metadata": "AppleDouble/macOS metadata, not Fortran source",
                        "test_driver": "test or example driver for a packaged surface",
                        "package_manifest": "package documentation stored with a Fortran suffix",
                        "helper_module": "support module without an independent PES entry",
                        "helper_source": "support source without an independent PES entry",
                        "adapter_module": "adapter module for another PES implementation",
                        "adapter": "coordinate/interface adapter for another PES implementation",
                        "coordinate_helper": "coordinate conversion helper",
                    }.get(classification, "non-surface package source"),
                    "owner_surface_ids": [],
                }
            )
        if read_units and not external:
            ambiguities.append(
                {
                    "path": record["path"],
                    "code": "unresolved_read_unit",
                    "detail": f"READ uses units {read_units} without a literal OPEN filename.",
                }
            )
        records[relative] = record

    identifiers: dict[str, list[str]] = defaultdict(list)
    for relative, record in records.items():
        if record["_surface"]:
            identifiers[record["id"]].append(relative)
    for duplicate_paths in identifiers.values():
        if len(duplicate_paths) < 2:
            continue
        for relative in duplicate_paths:
            parent = Path(relative).parent.name.lower()
            records[relative]["id"] = re.sub(
                r"[^a-z0-9]+", "-", parent
            ).strip("-")
    final_ids = [
        record["id"] for record in records.values() if record["_surface"]
    ]
    if len(final_ids) != len(set(final_ids)):
        raise ValueError("surface identifiers remain ambiguous after path disambiguation")

    build_dependencies(records)
    surface_ids = {
        relative: record["id"]
        for relative, record in records.items()
        if record["_surface"]
    }
    for relative, record in records.items():
        if record["_surface"]:
            continue
        record["owner_surface_ids"] = owner_ids(relative, records, surface_ids)
        if not record["owner_surface_ids"]:
            ambiguities.append(
                {
                    "path": record["path"],
                    "code": "missing_owner_surface",
                    "detail": "The archive contains no implementation source for this package manifest.",
                }
            )
    for record in records.values():
        for dependency in record["source_dependencies"]:
            dependency["path"] = f"src/{dependency['path']}"

    public_keys = {
        "_relative",
        "_active",
        "_declarations",
        "_surface",
    }
    surfaces = []
    non_surfaces = []
    for record in records.values():
        public = {key: value for key, value in record.items() if key not in public_keys}
        (surfaces if record["_surface"] else non_surfaces).append(public)

    classified_paths = [item["path"] for item in surfaces + non_surfaces]
    actual_surface_rejections = []
    for relative, record in records.items():
        if record["_surface"]:
            continue
        names = {
            item["name"].lower()
            for item in record["_declarations"]
            if item["kind"] in {"subroutine", "entry"}
        }
        fields = header_fields(read_source(src / relative))
        if names & {"pot", "surf"} or fields.get("interface", "").upper().startswith(
            ("HO-MM", "HE-MM")
        ):
            actual_surface_rejections.append(record["path"])
    validation = {
        "expected_fortran_files": EXPECTED_FILE_COUNT,
        "discovered_fortran_files": len(paths),
        "classified_fortran_files": len(classified_paths),
        "unique_classified_paths": len(set(classified_paths)),
        "all_files_classified_once": (
            len(paths) == len(classified_paths) == len(set(classified_paths))
        ),
        "matches_expected_vendor_count": len(paths) == EXPECTED_FILE_COUNT,
        "actual_surface_rejections": actual_surface_rejections,
        "no_actual_surface_rejected": not actual_surface_rejections,
        "unowned_non_surface_files": sum(
            not item["owner_surface_ids"] for item in non_surfaces
        ),
    }
    counts = {
        "fortran_files": len(paths),
        "primary_surface_sources": len(surfaces),
        "non_surface_files": len(non_surfaces),
        "ambiguities": len(ambiguities),
        "interface_families": dict(
            sorted(
                (
                    family,
                    sum(item["interface_family"] == family for item in surfaces),
                )
                for family in {item["interface_family"] for item in surfaces}
            )
        ),
    }
    return {
        "schema_version": 2,
        "vendor_root": "upstream/potlib2023",
        "counts": counts,
        "validation": validation,
        "ambiguities": sorted(
            ambiguities, key=lambda item: (item["path"], item["code"])
        ),
        "surfaces": surfaces,
        "non_surfaces": non_surfaces,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Inventory PotLib 2023 Fortran surface and support sources."
    )
    parser.add_argument("vendor_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--meson", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = inventory(args.vendor_root)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if args.meson:
        for surface in result["surfaces"]:
            atom_count = surface["atom_count"]
            test_atoms = 0 if atom_count["kind"] == "fixed" else 3
            print(f"{surface['id']}|{surface['path']}|{test_atoms}")
    validation = result["validation"]
    print(
        f"{result['counts']['fortran_files']} files: "
        f"{result['counts']['primary_surface_sources']} surfaces, "
        f"{result['counts']['non_surface_files']} non-surfaces, "
        f"{result['counts']['ambiguities']} ambiguities",
        file=sys.stderr,
    )
    if not (
        validation["all_files_classified_once"]
        and validation["matches_expected_vendor_count"]
        and validation["no_actual_surface_rejected"]
    ):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
