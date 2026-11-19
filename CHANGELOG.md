<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->
# Changelog

All notable changes to potlibc are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
  `fh2o-pipnn`, `soc-fh2o-nn`, `ho3-pes`, `h3o2-pipnn`, `pes-h3cl`,
  `ch3nh2`, and `pes-shell`.
- Conversion uses shared constants from `tools/generate_surface.py`:
  Bohr coordinates are multiplied by `BOHR_TO_ANGSTROM` before the
  deposited routine runs, kcal/mol energies are divided by
  `KCAL_PER_HARTREE`, and eV energies are divided by `EV_PER_HARTREE`. `ch3nh2` and `pes-shell` apply only the coordinate
  conversion because their deposited routines already return Hartree.

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
