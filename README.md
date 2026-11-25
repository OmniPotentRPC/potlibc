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

`test-tsan` is defined only when the C compiler accepts `-fsanitize=thread`
at compile time and at link time.

The runnable C examples are in [`examples/`](examples/); Meson builds and
executes them as part of the test suite.

Use `potlibc_create()` for fixed-size surfaces or `potlibc_create_sized()` for
variable-body surfaces with an id from [SURFACES.md](SURFACES.md), inspect
capabilities with `potlibc_get_metadata()`, and then call `potlibc_evaluate()`
or `potlibc_hessian_vector_product()`. HVP is generated with forward-dual
directional propagation, never a materialized Hessian, for every surface.
CBE (`ch52008`) also reports `POTLIBC_CAP_GRADIENT_JET`, and
`potlibc_gradient_jet()` returns the first three directional derivatives of
its analytic gradient with an operation-count roundoff report. The jet orders
are `H v`, `D^3 V(v, v)` and `D^4 V(v, v, v)`, laid out as
`[order][count][3 * atom_count]`. They come from the deposited analytic
gradient differentiated in forward mode, with no differencing step, and the
report holds the jet operation count and the largest magnitude among the
gradient and the returned orders. Their product with `DBL_EPSILON` is an
operation-count estimate of the roundoff, not a rigorous bound. Other surfaces
return `POTLIBC_ERROR_UNSUPPORTED` for it.
The all-surface consistency test uses
`potlibc_energy_directional_derivative()` as the dual oracle for each
analytic force, and a Richardson-extrapolated centered difference of that
derivative as the check on the projected Hessian-vector product. A nonzero
metadata reference geometry is only the probe center. CBE deposits one.
Every other surface reports a zero reference geometry, and the test
synthesizes a probe geometry.

The vendor tree retains the PotLib 2023 deposit apart from documented local
numerical corrections. `tools/surface_inventory.py` classifies the deposit and
`tools/generate_surface.py` derives family adapters from POTLIB headers and
interface metadata at build time. [SURFACES.md](SURFACES.md) explicitly accounts
for every Fortran source entry. [PROVENANCE.md](PROVENANCE.md) records every
vendor-file hash and correction.
The generated-source transforms and state-isolation design are documented in
[docs/design.md](docs/design.md).

The per-handle loader uses Linux `memfd_create` and `/proc/self/fd`.

See [`CHANGELOG.md`](CHANGELOG.md) for notable changes. To cite potlibc, use
the metadata in [`CITATION.cff`](CITATION.cff), and also cite the PotLib 2023
deposit that it lists under `references`.
