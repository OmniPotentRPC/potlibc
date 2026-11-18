/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#ifndef POTLIBC_PLUGIN_ABI_H
#define POTLIBC_PLUGIN_ABI_H

/* Plugins built before potlib_plugin_abi_version existed. */
#define POTLIBC_PLUGIN_ABI_LEGACY 1
/* Current generated-plugin ABI. Keep generate_surface.py in sync. */
#define POTLIBC_PLUGIN_ABI_VERSION 2
/* potlib_plugin_gradient_scale is optional; absence means this scale. */
#define POTLIBC_DEFAULT_GRADIENT_SCALE 1.0

/* Exported by current plugins and checked by libpotlibc at load. */
#define POTLIBC_PLUGIN_ABI_SYMBOL "potlib_plugin_abi_version"

#endif
