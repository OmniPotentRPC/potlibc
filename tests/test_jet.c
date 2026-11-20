/* SPDX-FileCopyrightText: 2026 Rohit Goswami */
/* SPDX-License-Identifier: BSD-3-Clause */

/* Bit goldens, report and argument checks of potlibc_gradient_jet on the CH5
 * surface. The goldens are the values of the standalone CBE jet kernel of the
 * spectral-instanton reference build: the jet orders, the operation count and
 * the scale of the roundoff report are all reproduced exactly. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "potlibc.h"
#include "ch52008_jet_oracle.h"

#define N 18

static int has_jet(potlibc_handle *handle)
{
    potlibc_metadata metadata;
    return potlibc_get_metadata(handle, &metadata) == POTLIBC_SUCCESS &&
           (metadata.capabilities & POTLIBC_CAP_GRADIENT_JET) != 0;
}

static int same_bits(const double *actual, const double *expected, size_t n)
{
    return memcmp(actual, expected, n * sizeof(double)) == 0;
}

static int check_jet_at(potlibc_handle *handle, const double *x,
                        const double *direction, double energy,
                        const double *gradient, const double *hvp,
                        const double *cubic, const double *quartic,
                        long operations, double scale, long hvp_operations,
                        double hvp_scale)
{
    double e;
    double force[N];
    double jet[3 * 2 * N];
    double doubled_x[2 * N];
    double doubled_direction[2 * N];
    double product[N];
    potlibc_jet_report report = {0, 0.0};
    size_t i;

    if (potlibc_evaluate(handle, 1, x, &e, force) != POTLIBC_SUCCESS ||
        memcmp(&e, &energy, sizeof e) != 0) {
        return 1;
    }
    for (i = 0; i < N; i++) {
        const double negated = -gradient[i];
        if (memcmp(&negated, &force[i], sizeof negated) != 0) {
            return 2;
        }
    }
    if (potlibc_gradient_jet(handle, 1, 3, x, direction, jet, &report) !=
        POTLIBC_SUCCESS) {
        return 3;
    }
    if (!same_bits(jet, hvp, N) || !same_bits(jet + N, cubic, N) ||
        !same_bits(jet + 2 * N, quartic, N)) {
        return 4;
    }
    if (report.operations != operations ||
        memcmp(&report.scale, &scale, sizeof scale) != 0) {
        return 5;
    }
    if (potlibc_gradient_jet(handle, 1, 1, x, direction, jet, &report) !=
            POTLIBC_SUCCESS ||
        !same_bits(jet, hvp, N) || report.operations != hvp_operations ||
        memcmp(&report.scale, &hvp_scale, sizeof hvp_scale) != 0) {
        return 6;
    }
    if (potlibc_gradient_jet(handle, 1, 3, x, direction, jet, NULL) !=
        POTLIBC_SUCCESS) {
        return 7;
    }
    memcpy(doubled_x, x, sizeof(double) * N);
    memcpy(doubled_x + N, x, sizeof(double) * N);
    memcpy(doubled_direction, direction, sizeof(double) * N);
    memcpy(doubled_direction + N, direction, sizeof(double) * N);
    if (potlibc_gradient_jet(handle, 2, 3, doubled_x, doubled_direction, jet,
                             &report) != POTLIBC_SUCCESS ||
        report.operations != 2 * operations ||
        !same_bits(jet, hvp, N) || !same_bits(jet + N, hvp, N) ||
        !same_bits(jet + 2 * N, cubic, N) ||
        !same_bits(jet + 3 * N, cubic, N) ||
        !same_bits(jet + 4 * N, quartic, N) ||
        !same_bits(jet + 5 * N, quartic, N)) {
        return 8;
    }
    if (potlibc_hessian_vector_product(handle, 1, x, direction, product) !=
            POTLIBC_SUCCESS ||
        !same_bits(product, hvp, N)) {
        return 9;
    }
    return 0;
}

int main(void)
{
    potlibc_handle *handle = NULL;
    potlibc_handle *proton = NULL;
    double jet[3 * N];
    double bad[N];
    potlibc_jet_report report = {-1, -1.0};
    int failure;

    if (potlibc_create("ch52008", &handle) != POTLIBC_SUCCESS) return 1;
    if (!has_jet(handle)) {
        return 2;
    }
    if (potlibc_gradient_jet(handle, 1, 0, saddle_x, saddle_v, jet, &report) !=
            POTLIBC_ERROR_JET_ORDER ||
        potlibc_gradient_jet(handle, 1, 4, saddle_x, saddle_v, jet, &report) !=
            POTLIBC_ERROR_JET_ORDER ||
        report.operations != -1) {
        return 3;
    }
    if (potlibc_gradient_jet(handle, 0, 1, saddle_x, saddle_v, jet, &report) !=
        POTLIBC_ERROR_INVALID_COUNT) {
        return 4;
    }
    if (potlibc_gradient_jet(handle, 1, 1, NULL, saddle_v, jet, &report) !=
            POTLIBC_ERROR_NULL_POINTER ||
        potlibc_gradient_jet(handle, 1, 1, saddle_x, NULL, jet, &report) !=
            POTLIBC_ERROR_NULL_POINTER ||
        potlibc_gradient_jet(handle, 1, 1, saddle_x, saddle_v, NULL,
                             &report) != POTLIBC_ERROR_NULL_POINTER ||
        potlibc_gradient_jet(NULL, 1, 1, saddle_x, saddle_v, jet, &report) !=
            POTLIBC_ERROR_NULL_POINTER) {
        return 5;
    }
    memcpy(bad, saddle_v, sizeof bad);
    bad[3] = NAN;
    if (potlibc_gradient_jet(handle, 1, 1, saddle_x, bad, jet, &report) !=
        POTLIBC_ERROR_NONFINITE_INPUT) {
        return 6;
    }
    failure = check_jet_at(handle, saddle_x, saddle_v, saddle_energy,
                           saddle_gradient, saddle_hvp, saddle_cubic,
                           saddle_quartic, saddle_jet_operations,
                           saddle_jet_scale, saddle_hvp_operations,
                           saddle_hvp_scale);
    if (failure != 0) {
        fprintf(stderr, "saddle jet golden failed: %d\n", failure);
        return 20 + failure;
    }
    failure = check_jet_at(handle, displaced_x, displaced_v, displaced_energy,
                           displaced_gradient, displaced_hvp, displaced_cubic,
                           displaced_quartic, displaced_jet_operations,
                           displaced_jet_scale, displaced_hvp_operations,
                           displaced_hvp_scale);
    if (failure != 0) {
        fprintf(stderr, "displaced jet golden failed: %d\n", failure);
        return 40 + failure;
    }
    potlibc_destroy(handle);

    /* A surface without the jet kernel reports it as unsupported. */
    if (potlibc_create("chc", &proton) == POTLIBC_SUCCESS) {
        if (has_jet(proton) ||
            potlibc_gradient_jet(proton, 1, 1, saddle_x, saddle_v, jet,
                                 &report) != POTLIBC_ERROR_UNSUPPORTED) {
            return 60;
        }
        potlibc_destroy(proton);
    }
    return 0;
}
