<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->
# Generated sources and state isolation

`upstream/potlib2023/` retains the PotLib 2023 deposit apart from documented
local numerical corrections. The build compiles the corrected sources through
derived files written into the Meson build directory.
`tools/surface_inventory.py` classifies all 403 deposited Fortran files.
`tools/generate_surface.py` writes a generated surface, utility unit, and C ABI
adapter into the Meson build directory for each of the 387 surface roots.

The source transform normalizes legacy fixed-form declarations, repairs
deposited declarations that modern gfortran rejects, embeds parameter files
into per-plugin initialization, and generates the adapter selected by the
declared interface family or callable signature. Mixed forward-dual values
carry a coordinate seed, a direction seed, and their mixed term. Forces read
the coordinate derivative and HVPs read the mixed derivative, so no full
Hessian is generated or stored. CBE retains the deposited analytic force path
for bit parity and differentiates that analytic gradient in the direction
seed for its HVP.

The CBE surface adds a third generated variant on third-order jets
(`src/jet3_mod.f90`): a value and its first three derivatives along one
direction, with exact product rules. `tools/catalog.py::to_jet3` retypes the
`dual2` variant under `j`-prefixed names, and wraps every double-precision
literal outside a `PARAMETER` as a jet constant so that a constant operand is
counted like any other jet operation. The adapter exports the optional symbol
`potlib_plugin_gradient_jet`, which `potlibc_gradient_jet()` looks up to set
`POTLIBC_CAP_GRADIENT_JET`. Only `ch52008` carries it, so no other plugin links
the jet module. The jet is compiled with the surface at `-O3` and reproduces
the standalone CBE kernel bit for bit, including its operation count.

The dual intrinsics use explicit boundary conventions. `abs(0)` and `sqrt(0)`
select a zero subgradient, constant `exp` tangents remain zero even when the
primal overflows, and integral-valued real powers use the integer rule so
negative bases do not pass through `log`. Mixed scalar/dual `min` and `max`
overloads preserve a constant bound at a clipping tie. Comparisons otherwise
select branches from the primal value; consistency probes avoid
nondifferentiable switching surfaces.

`ch52008` is compiled at `-O3`. `phsch3-ia`, `phsch3-s-1-8`,
`phsch3-s-2-0`, `phsch3-ii`, `phsch3-s-2-5`, `phsch3-iii`,
`phsch3-s-2-2`, `phsch3-aprp`, and `ch4t8` are compiled at `-O0`, because
higher optimization faults those legacy implementations even though the
generated adapter is unchanged. `-finit-derived` and `-finit-real=zero`
are applied to `h3dmbe`, `sioh`, and `oh2m3a2001`; to that `phsch3` set
plus `ch3nh2`, `phsh-aprp`, and `ch4t8`; and to every surface id that
starts with `c2h6`. These are per-surface compiler rules, not edits to
the vendor sources.

POTLIB surfaces retain mutable process state in COMMON blocks and variables
with SAVE storage. Every generated plugin is linked with symbolic references
bound to that plugin. `potlibc_create` copies the plugin bytes into a new
anonymous file and loads that file with `dlopen`. The dynamic linker therefore
creates a distinct data image for every handle, while shared dependencies such
as libgfortran and libc remain in the base namespace. Closing a handle unloads
only its private plugin image.

After initialization and metadata selection, every generated adapter other
than CBE (`ch52008`) snapshots its writable ELF data segments. Each
coordinate or direction seed restores that snapshot before entering the
transformed surface. CBE keeps the deposited analytic `pot` path and does
not restore a snapshot. Its isolation is still the private plugin image.
The snapshot preserves initialized parameter tables while preventing work
arrays, COMMON blocks, and SAVE variables mutated by one seed from
contaminating another seed. Formatted parameter records are read into the
real lane of dual values and their derivative lanes are initialized
explicitly.

This design uses no global lock and does not consume glibc ELF namespaces.
glibc's thread-safe `dlopen` and `dlclose` implementation handles concurrent
image creation and removal; evaluations touch only their handle's private
image and, except for CBE, its snapshot. The parallel test keeps 128 handles
live across 64 threads and compares six surfaces bit-for-bit with serial
evaluation. That test is also built and run with ThreadSanitizer when the
compiler accepts `-fsanitize=thread` at compile time and at link time.
