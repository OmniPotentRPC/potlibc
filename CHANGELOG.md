<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->
# Changelog

All notable changes to potlibc are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.2.1] - 2026-11-21

### Fixed

- Every math-using C test now links `libm`, which `tests/meson.build` looks up
  once and keeps optional. The symmetric `dual2` Hessian test did not link it
  at all, so it could fail to link on toolchains where `libm` is a separate
  library.
- The docs linkcheck ignores `iopscience.iop.org` and sets
  `linkcheck_timeout = 20` and `linkcheck_retries = 2`, so a slow host or a
  transient network error no longer fails the build.
- The documented total in `docs/source/accuracy.rst` is 403 of 403 Meson tests
  (previously 400), with 387 of 387 per-surface dual consistency cases.
- `tools/validate_generated_surfaces.py` compiles `dual2_mod` and `jet3_mod`
  and links each generated surface plugin to the modules that its `use`
  statements name, instead of always linking `dual2_mod`. `ch52008` uses
  `jet3_mod` in its generated source, utility, and adapter, so it did not
  compile before. This change shipped in 0.2.0 but was missing from the 0.2.0
  notes.

## [0.2.0] - 2026-11-21

### Added

- `potlibc_gradient_jet()` and `POTLIBC_CAP_GRADIENT_JET`: the first three
  directional derivatives of the CBE (`ch52008`) analytic gradient, laid out
  as `[order][count][3 * atom_count]`, with a `potlibc_jet_report` of the jet
  operation count and the scale for an operation-count roundoff estimate.
  Other surfaces return `POTLIBC_ERROR_UNSUPPORTED`. The new status
  `POTLIBC_ERROR_JET_ORDER` is appended to the enum.
- `jet3_mod`, a third-order forward-mode jet, and a generated `jet3` variant of
  the CH5 surface. Tests pin its bits, operation count and scale against the
  standalone CBE kernel, check each order against a Richardson-extrapolated
  finite difference of the lower orders, require that a 1e-5 relative error in
  an order-3 component is detected, and pin integer powers at zero.

### Fixed

- The files that `docs/generate_surface_catalog.py` and `tools/catalog.py`
  generate carry valid SPDX headers. The copyright text and the license
  expression each sit on their own line, so the license reads
  `BSD-3-Clause` instead of including Python string quotes.

## [0.1.1] - 2026-11-20

### Fixed

- Preserved the LSTH and ALC auxiliary-barrier Hessian on their symmetric
  coordinate planes by squaring the signed permutation factor before applying
  its absolute-value convention to odd powers.
- Mixed scalar and dual `min` and `max` keep a constant bound, and a zero
  derivative, when the primal values tie. This includes a dual argument
  followed by a real or integer bound.
- The following adapters now return native Bohr and Hartree values instead
  of the v0.1.0 behavior: `n3pes`, `o3pes`, `o4pes`, `c2h7`,
  `c2h7pot5`, `no2-2ap-pip-nn`, `no2-4ap-pip-nn`, `no2-6ap-pip-nn`,
  `fh2o-pipnn`, `soc-fh2o-nn`, `ho3-pes`, `h3o2-pipnn`, `pes-h3cl`, and
  `ch3nh2`.
- Conversion uses shared constants from `tools/generate_surface.py`:
  Bohr coordinates are multiplied by `BOHR_TO_ANGSTROM` before the
  deposited routine runs, kcal/mol energies are divided by
  `KCAL_PER_HARTREE`, and eV energies are divided by `EV_PER_HARTREE`.
  `ch3nh2` applies only the coordinate conversion because its deposited
  routine already returns Hartree.
- `pes-shell` forwards Bohr coordinates directly to `f` instead of
  multiplying by `BOHR_TO_ANGSTROM`, because `pes_shell.f90` uses
  `a0 = 2.d0` Bohr (`test.f90` divides input Angstrom coordinates by
  `auang` before calling `f`) and already returns Hartree energies and
  Hartree/Bohr gradients.
- The public header, README, and user guide now match the implementation.
  `POTLIBC_ERROR_INVALID_SURFACE_ID` is a syntax and path check. A well
  formed identifier with no loadable plugin returns
  `POTLIBC_ERROR_SURFACE_UNAVAILABLE`. `potlibc_create_sized` rejects a
  zero atom count with `POTLIBC_ERROR_INVALID_ATOM_COUNT`. The consistency
  test uses `potlibc_energy_directional_derivative()` as the force oracle.
  CBE (`ch52008`) does not restore a writable data snapshot. Documented
  compiler flags match `meson.build`, and `test-tsan` exists only when the
  compiler accepts `-fsanitize=thread` at compile time and at link time.
  The regression suite count is 400.

## [0.1.0] - 2026-11-18

### Added

- Vendored the PotLib 2023 surface deposit verbatim, with a recorded
  SHA-256 digest for every vendor file.
- A stable C ABI (`potlibc_create`, `potlibc_create_sized`,
  `potlibc_evaluate`, `potlibc_energy_directional_derivative`,
  `potlibc_hessian_vector_product`, `potlibc_destroy`) and the Meson build
  that generates one isolated plugin per catalog surface.
- Generated per-surface Fortran adapters, parameter data, and a
  mixed-forward-dual implementation so every surface exposes energy,
  analytic forces, and an analytic Hessian-vector product without forming
  a dense Hessian.
- A parity and parallel-isolation test suite, runnable C documentation
  examples, and a Sphinx and Doxygen user and API guide covering
  provenance, generated-source design, accuracy, and the full surface
  catalog.
- A documentation link check and a DOI resolution check, run alongside the
  existing warnings-as-errors documentation build and regression suite.

### Changed

- Enriched surface catalog metadata and recognized additional vendor
  header variants when extracting each surface's formula, coordinates,
  reference, and DOIs.

### Fixed

- Linked the C API test against `libm` so Meson builds it without a
  missing symbol on toolchains that do not link it implicitly.
- Derived each surface's internal-coordinate ordering and metadata
  directly from its declared interface family instead of a shared
  default, and honored the per-surface triatomic internal-coordinate
  order this requires.
- Replaced the GitHub Pages workflow's unresolvable pinned action
  references, which had been silently failing every deployment to
  `main`, with the actions' real release commits.
- Corrected the test count documented for the regression suite, which had
  fallen out of date after two triatomic internal-coordinate tests were
  added.
