.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

rgpot integration status
========================

`rgpot <https://rgpot.rgoswami.me>`_ is the ecosystem landing page for
potential selection, simulation-engine integrations, and language bindings.
The direct ``potlibc`` C ABI documented here is complete, but rgpot commit
``672aef430855036125c902ace024756d682e180a`` does **not** yet contain a
``potlibc``-backed optional potential.

Consequently, there is currently no truthful rgpot C++ construction example,
Python class, or Rust wrapper to document:

* rgpot's C++ potential hierarchy exists, but has no ``potlibc`` adapter.
* rgpot's Python nanobind package exposes selected rgpot potentials, but no
  ``potlibc`` surface selector.
* ``rgpot-core`` is a Rust crate, but it does not bind ``potlibc``.

Use the C ABI directly until those integrations land. The following issue text
records the gaps without reserving API names.

TODO issue: optional rgpot C++ potential
----------------------------------------

**Title:** Add an optional rgpot potential backed by potlibc

**Scope:** Add a Meson feature option that discovers ``potlibc`` through
pkg-config, implement an rgpot C++ ``Potential`` adapter with explicit unit
conversion and surface/atom-count configuration, and register it only when
the dependency is enabled.

**Acceptance:** Energy and force batches agree with direct C calls; HVP
capability is forwarded; one independent ``potlibc_handle`` is used per
concurrent evaluator; disabled builds retain the same ABI and dependency set.

TODO issue: expose the adapter in Python
----------------------------------------

**Title:** Bind the optional rgpot potlibc potential in Python

**Blocked by:** the C++ adapter issue above.

**Acceptance:** The generated stubs, wheel-content tests, and Python API
reference cover only symbols built with the optional dependency; unavailable
builds fail with an actionable capability message.

TODO issue: decide Rust binding ownership
-----------------------------------------

**Title:** Define whether potlibc belongs in rgpot-core

**Scope:** Decide between a small dedicated ``potlibc-sys`` crate and an
optional rgpot-core feature before exposing a safe wrapper. Specify handle
ownership, ``Send``/``Sync`` behavior, plugin discovery, and version checks.

**Acceptance:** No safe concurrency traits are declared without tests proving
their contracts, and the wrapper is tested against the installed C ABI.
