#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
"""Generate, compile, link, and inspect every inventoried surface plugin."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True

try:
    from .generate_surface import _adapter, _plan, _transform
    from .surface_inventory import inventory
except ImportError:
    from generate_surface import _adapter, _plan, _transform
    from surface_inventory import inventory


REQUIRED_SYMBOLS = {
    "potlib_plugin_initialize",
    "potlib_plugin_atom_count",
    "potlib_plugin_metadata",
    "potlib_plugin_evaluate",
    "potlib_plugin_energy_directional",
    "potlib_plugin_hvp",
    "potlibc_state_checkpoint",
    "potlibc_state_restore",
    "potlibc_state_release",
}
FORTRAN_FLAGS = (
    "-shared",
    "-fPIC",
    "-std=legacy",
    "-ffree-line-length-none",
    "-fallow-argument-mismatch",
    "-fno-align-commons",
)
FORTRAN_MODULES = {
    "dual2_mod": "dual2_mod.f90",
    "jet3_mod": "jet3_mod.f90",
}
_USE_MODULE = re.compile(
    r"(?im)^[ \t]*use[ \t]+([A-Za-z_][A-Za-z0-9_]*)\b"
)


def _tree_digest(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _modules_used(directory: Path) -> list[str]:
    names: set[str] = set()
    for filename in ("source.f90", "utility.f90", "adapter.f90"):
        text = (directory / filename).read_text()
        names.update(match.group(1).lower() for match in _USE_MODULE.finditer(text))
    return [name for name in FORTRAN_MODULES if name in names]


def _compile(item: dict[str, Any]) -> dict[str, Any]:
    directory = Path(item["directory"])
    command = [
        "gfortran",
        *FORTRAN_FLAGS,
        f"-I{item['module_directory']}",
        "-J.",
        "source.f90",
        "utility.f90",
        "adapter.f90",
        *item["module_objects"],
        item["state_object"],
        item["data_stub_object"],
        "-ldl",
        "-llapack",
        "-lblas",
        "-Wl,-z,defs",
        "-o",
        "plugin.so",
    ]
    result = subprocess.run(
        command, cwd=directory, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False
    )
    record: dict[str, Any] = {
        "id": item["id"],
        "path": item["path"],
        "command": command,
        "compiled": result.returncode == 0,
    }
    if result.returncode != 0:
        record["cause"] = result.stdout.strip()
        return record
    symbols = {
        line.split()[-1]
        for line in subprocess.check_output(
            ["nm", "-D", "--defined-only", "plugin.so"],
            cwd=directory,
            text=True,
        ).splitlines()
        if line.split()
    }
    missing = sorted(REQUIRED_SYMBOLS - symbols)
    record["linked"] = True
    record["symbols"] = not missing
    if missing:
        record["cause"] = f"missing ABI symbols: {', '.join(missing)}"
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("vendor", type=Path)
    parser.add_argument("build_directory", type=Path)
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()
    vendor = args.vendor.resolve()
    build = args.build_directory.resolve()
    build.mkdir(parents=True, exist_ok=True)
    module_directory = build / "modules"
    module_directory.mkdir(exist_ok=True)
    before = _tree_digest(vendor)
    source_root = Path(__file__).resolve().parents[1] / "src"
    module_objects = {}
    module_commands = []
    for name, filename in FORTRAN_MODULES.items():
        object_path = build / f"{name}.o"
        command = [
            "gfortran", "-c", "-fPIC", f"-J{module_directory}",
            str(source_root / filename),
            "-o", str(object_path),
        ]
        subprocess.run(command, check=True)
        module_objects[name] = object_path
        module_commands.append(command)
    state_object = build / "plugin_state.o"
    state_command = [
        "cc", "-c", "-fPIC",
        str(Path(__file__).resolve().parents[1] / "src" / "plugin_state.c"),
        "-o", str(state_object),
    ]
    subprocess.run(state_command, check=True)
    data_stub_source = build / "data_stub.c"
    data_stub_source.write_text(
        "int potlibc_n4_data_fd(void) { return -1; }\n"
    )
    data_stub_object = build / "data_stub.o"
    data_stub_command = [
        "cc", "-c", "-fPIC", str(data_stub_source),
        "-o", str(data_stub_object),
    ]
    subprocess.run(data_stub_command, check=True)

    catalog = inventory(vendor)
    generated: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for index, record in enumerate(catalog["surfaces"], 1):
        directory = build / f"{index:03d}-{record['id']}"
        directory.mkdir(exist_ok=True)
        try:
            plan = _plan(record, vendor)
            source, utility = _transform(plan, vendor)
            (directory / "source.f90").write_text(source)
            (directory / "utility.f90").write_text(utility)
            (directory / "adapter.f90").write_text(_adapter(plan))
            used = _modules_used(directory)
        except Exception as error:
            failures.append(
                {
                    "id": record["id"],
                    "path": record["path"],
                    "generated": False,
                    "cause": f"{type(error).__name__}: {error}",
                }
            )
            continue
        generated.append(
            {
                "id": record["id"],
                "path": record["path"],
                "directory": str(directory),
                "module_directory": str(module_directory),
                "module_objects": [str(module_objects[name]) for name in used],
                "state_object": str(state_object),
                "data_stub_object": str(data_stub_object),
            }
        )

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(_compile, generated))
    failures.extend(item for item in results if not item.get("symbols"))
    after = _tree_digest(vendor)
    summary = {
        "inventoried": len(catalog["surfaces"]),
        "generated": len(generated),
        "compiled": sum(bool(item.get("compiled")) for item in results),
        "linked": sum(bool(item.get("linked")) for item in results),
        "abi_complete": sum(bool(item.get("symbols")) for item in results),
        "vendor_unchanged": before == after,
        "vendor_sha256": after,
        "module_commands": module_commands,
        "state_command": state_command,
        "data_stub_command": data_stub_command,
        "fortran_flags": list(FORTRAN_FLAGS),
        "failures": failures,
    }
    (build / "report.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items()
                      if key != "failures"}, indent=2))
    print(f"failures: {len(failures)}")
    return 0 if not failures and before == after else 1


if __name__ == "__main__":
    raise SystemExit(main())
