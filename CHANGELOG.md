<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->
# Changelog

All notable changes to potlibc are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
