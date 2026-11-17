<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->
# Generated sources and state isolation

`upstream/potlib2023/` is an immutable copy of the PotLib 2023 deposit. The
build never compiles modified files in that directory. Instead,
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

The dual intrinsics use explicit boundary conventions. `abs(0)` and `sqrt(0)`
select a zero subgradient, constant `exp` tangents remain zero even when the
primal overflows, and integral-valued real powers use the integer rule so
negative bases do not pass through `log`. Mixed scalar/dual `min` and `max`
overloads preserve a constant bound at a clipping tie. Comparisons otherwise
select branches from the primal value; consistency probes avoid
nondifferentiable switching surfaces.

The build applies zero initialization only to the small set of deposited
routines that leave transformed work values unset. It also compiles the
optimization-sensitive APRP and CH4T8 routines at `-O0`; higher optimization
causes those legacy implementations to fault even though the generated
adapter is unchanged. These are per-surface compiler rules, not edits to the
vendor sources.

POTLIB surfaces retain mutable process state in COMMON blocks and variables
with SAVE storage. Every generated plugin is linked with symbolic references
bound to that plugin. `potlibc_create` copies the plugin bytes into a new
anonymous file and loads that file with `dlopen`. The dynamic linker therefore
creates a distinct data image for every handle, while shared dependencies such
as libgfortran and libc remain in the base namespace. Closing a handle unloads
only its private plugin image.

After initialization and metadata selection, the plugin snapshots its writable
ELF data segments. Each coordinate or direction seed restores that snapshot
before entering the transformed surface. This preserves initialized parameter
tables while preventing work arrays, COMMON blocks, and SAVE variables mutated
by one seed from contaminating another seed. Formatted parameter records are
read into the real lane of dual values and their derivative lanes are
initialized explicitly.

This design uses no global lock and does not consume glibc ELF namespaces.
glibc's thread-safe `dlopen` and `dlclose` implementation handles concurrent
image creation and removal; evaluations touch only their handle's private
image and snapshot. The parallel test keeps 128 handles live across 64 threads
and compares six surfaces bit-for-bit with serial evaluation. The same test is
also built and run with ThreadSanitizer.
