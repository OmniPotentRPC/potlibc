<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->

# Changelog

## 0.2.0

Corrected CH5 2008 CBE forces and Hessian-vector products are the default.
`potlibc_set_gradient_mode()` selects `POTLIBC_GRADIENT_RAW_UPSTREAM` when an
application needs the deposited gradient bits. The public C ABI adds that call
and the `potlibc_gradient_mode` enum. Creation and evaluation signatures are
unchanged, and the libpotlibc soname stays `libpotlibc.so.0`.

Current plugins export `potlib_plugin_abi_version`. An older plugin that omits
it still loads, and a missing `potlib_plugin_gradient_scale` uses scale 1.0.
Calling that version symbol records host acceptance. Until the call,
`potlib_plugin_atom_count` returns a negative count, so an older library
reports a plugin ABI mismatch instead of returning an uncorrected CH5 2008
gradient.
