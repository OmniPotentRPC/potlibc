/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include "potlibc.h"
#include "../src/plugin_abi.h"

#include <stdio.h>

static const double stub_energy = 1.5;
static const double stub_force_x = 2.0;
static const double stub_force_y = -3.0;
static const double stub_force_z = 0.5;
static const double current_gradient_scale = 2.0;

static int forces_match(const double force[3], double scale)
{
    return force[0] == stub_force_x * scale &&
           force[1] == stub_force_y * scale &&
           force[2] == stub_force_z * scale;
}

static int load_with_scale(const char *surface_id, double scale)
{
    potlibc_handle *handle = NULL;
    const double coordinates[3] = {0.1, 0.2, 0.3};
    double energy = 0.0;
    double force[3] = {0.0, 0.0, 0.0};
    potlibc_status status = potlibc_create(surface_id, &handle);
    if (status != POTLIBC_SUCCESS || handle == NULL) {
        fprintf(stderr, "%s: create returned %s\n", surface_id,
                potlibc_status_string(status));
        potlibc_destroy(handle);
        return 0;
    }
    status = potlibc_evaluate(handle, 1, coordinates, &energy, force);
    if (status != POTLIBC_SUCCESS || energy != stub_energy ||
        !forces_match(force, scale)) {
        fprintf(stderr,
                "%s: corrected force mismatch status=%s energy=%a "
                "force=%a %a %a scale=%a\n",
                surface_id, potlibc_status_string(status), energy, force[0],
                force[1], force[2], scale);
        potlibc_destroy(handle);
        return 0;
    }
    if (scale != POTLIBC_DEFAULT_GRADIENT_SCALE) {
        double raw_force[3];
        status = potlibc_set_gradient_mode(handle,
                                           POTLIBC_GRADIENT_RAW_UPSTREAM);
        if (status != POTLIBC_SUCCESS) {
            fprintf(stderr, "%s: raw mode rejected\n", surface_id);
            potlibc_destroy(handle);
            return 0;
        }
        status = potlibc_evaluate(handle, 1, coordinates, &energy, raw_force);
        if (status != POTLIBC_SUCCESS || energy != stub_energy ||
            !forces_match(raw_force, POTLIBC_DEFAULT_GRADIENT_SCALE)) {
            fprintf(stderr, "%s: raw mode changed the deposited force\n",
                    surface_id);
            potlibc_destroy(handle);
            return 0;
        }
    }
    potlibc_destroy(handle);
    return 1;
}

static int reject_future_plugin(void)
{
    potlibc_handle *handle = NULL;
    potlibc_status status = potlibc_create("abi-future", &handle);
    if (status != POTLIBC_ERROR_PLUGIN_ABI || handle != NULL) {
        fprintf(stderr, "abi-future: expected plugin ABI mismatch, got %s\n",
                potlibc_status_string(status));
        potlibc_destroy(handle);
        return 0;
    }
    return 1;
}

int main(void)
{
    if (!load_with_scale("abi-old", POTLIBC_DEFAULT_GRADIENT_SCALE)) {
        return 1;
    }
    if (!load_with_scale("abi-current", current_gradient_scale)) {
        return 2;
    }
    if (!reject_future_plugin()) {
        return 3;
    }
    return 0;
}
