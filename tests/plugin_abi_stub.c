/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include "../src/plugin_abi.h"

#if defined(PLUGIN_ABI_OLD)
#define PLUGIN_KIND_OLD 1
#elif defined(PLUGIN_ABI_CURRENT)
#define PLUGIN_KIND_CURRENT 1
#elif defined(PLUGIN_ABI_FUTURE)
#define PLUGIN_KIND_FUTURE 1
#else
#error "define PLUGIN_ABI_OLD, PLUGIN_ABI_CURRENT, or PLUGIN_ABI_FUTURE"
#endif

#if defined(PLUGIN_KIND_CURRENT)
static const double stub_gradient_scale = 2.0;
#endif

#if !defined(PLUGIN_KIND_OLD)
static int host_abi_accepted;
#endif

#if defined(PLUGIN_KIND_CURRENT)
int potlib_plugin_abi_version(void)
{
    host_abi_accepted = 1;
    return POTLIBC_PLUGIN_ABI_VERSION;
}

double potlib_plugin_gradient_scale(void)
{
    return stub_gradient_scale;
}
#elif defined(PLUGIN_KIND_FUTURE)
int potlib_plugin_abi_version(void)
{
    host_abi_accepted = 1;
    return POTLIBC_PLUGIN_ABI_VERSION + 1;
}
#endif

void potlib_plugin_initialize(void) {}

int potlib_plugin_atom_count(void)
{
#if defined(PLUGIN_KIND_OLD)
    return 1;
#else
    if (!host_abi_accepted) {
        return -1;
    }
    return 1;
#endif
}

int potlib_plugin_metadata(int atom_count, int *derivatives,
                           int *atomic_numbers, double *reference)
{
    if (atom_count != 1) {
        return 1;
    }
    *derivatives = 1;
    atomic_numbers[0] = 1;
    reference[0] = 0.0;
    reference[1] = 0.0;
    reference[2] = 0.0;
    return 0;
}

void potlib_plugin_evaluate(int count, const double *x, double *energy,
                            double *force)
{
    (void)x;
    for (int sample = 0; sample < count; ++sample) {
        energy[sample] = 1.5;
        force[3 * sample] = 2.0;
        force[3 * sample + 1] = -3.0;
        force[3 * sample + 2] = 0.5;
    }
}
