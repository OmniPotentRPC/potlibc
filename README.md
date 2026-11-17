<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->
# potlibc

potlibc provides a stable C ABI for all 387 surface versions in the PotLib
2023 deposit. Each handle loads a private in-memory
copy of its generated surface plugin, isolating Fortran COMMON and SAVE state
without a global lock or an ELF namespace per handle. Coordinates are Bohr,
energies Hartree, and forces Hartree/Bohr.

The full user guide, generated C API reference, and exhaustive surface catalog
are published at <https://omnipotentrpc.github.io/potlibc/>. The current
per-handle loader requires Linux.

```sh
meson setup build
meson compile -C build
meson test -C build
meson compile -C build test-tsan
```

The runnable C examples are in [`examples/`](examples/); Meson builds and
executes them as part of the test suite.

Use `potlibc_create()` for fixed-size surfaces or `potlibc_create_sized()` for
variable-body surfaces with an id from [SURFACES.md](SURFACES.md), inspect
capabilities with `potlibc_get_metadata()`, and then call `potlibc_evaluate()`
or `potlibc_hessian_vector_product()`. HVP is generated with forward-dual
directional propagation, never a materialized Hessian, for every surface.
The metadata's reference geometry and
`potlibc_energy_directional_derivative()` provide the independent dual oracle
used by the all-surface consistency test.

The vendor tree remains byte-for-byte unchanged. `tools/surface_inventory.py`
classifies the deposit and `tools/generate_surface.py` derives family adapters
from POTLIB headers and interface metadata at build time. [SURFACES.md](SURFACES.md)
explicitly accounts for every Fortran source entry. [PROVENANCE.md](PROVENANCE.md)
records every vendor-file hash.
The generated-source transforms and state-isolation design are documented in
[docs/design.md](docs/design.md).

The per-handle loader uses Linux `memfd_create` and `/proc/self/fd`.
