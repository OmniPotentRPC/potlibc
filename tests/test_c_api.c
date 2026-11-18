/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include "potlibc.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static const double CBE_GRADIENT_SCALE = 0x1.000022993a1a7p+0;

static const double x[18] = {
    2.5448986937601746, 3.0319878740728727, 2.603957609783651,
    1.1100980488026544, 1.3613545004038685, 1.1486133528262976,
    0.013880929393421192, 0.2086334261632169, 2.457516169583444,
    2.4611712350631576, 0.21486419033998002, 0.07550700261319757,
    -0.09665148160819977, 2.4937031206965017, -0.0782634154328256,
    3.3853544918540113, 4.100244654765713, 3.5175358525837166,
};
static const double direction[18] = {
    -0.31, 0.29, -0.23, 0.19, -0.17, 0.13, -0.11, 0.07, -0.05,
    0.03, -0.02, 0.01, 0.04, -0.06, 0.08, -0.12, 0.16, -0.18,
};
static const double expected_energy = 0x1.885b2451b0c60p-6;
static const double expected_force[18] = {
    -0x1.8316b6e1e082bp-12, 0x1.1201ad10d0427p-10,
    0x1.2df60c85d1378p-11, 0x1.28b85073adbadp-9,
    -0x1.04481763e2e94p-8, -0x1.4069c0e61ceccp-9,
    -0x1.4cc897c6a6f6dp-13, 0x1.8312a8eb50d7cp-12,
    0x1.977fd446af83dp-12, -0x1.49438572248c4p-9,
    0x1.41bb2fdbe7984p-9, 0x1.7defd11a38807p-10,
    0x1.9220ac0c86087p-12, 0x1.9eec4f2507d27p-13,
    0x1.a0e04b091d913p-19, 0x1.9bb3feac645efp-12,
    -0x1.8f9e358890dfcp-14, 0x1.4e114dfa1455fp-16,
};
static const double expected_hvp[18] = {
    -0x1.7ddbef4065ddcp-11, 0x1.997d00896570cp-6,
    0x1.95c8d5837a79cp-9, 0x1.1e7f15301e18dp-4,
    -0x1.9c2884f59d46fp-4, 0x1.67ba3cb02433dp-5,
    -0x1.8c6d9a9cfce27p-10, 0x1.b8cfc35517323p-6,
    -0x1.86658c96a6ee4p-6, -0x1.bdfab35359e1ap-6,
    0x1.c4470493f4e69p-6, 0x1.215ed7e63d170p-7,
    -0x1.2aa7d2c4dfedfp-5, 0x1.09a47ff7e4a43p-5,
    -0x1.aaec0d1de68f5p-6, -0x1.07f10afe2ea37p-8,
    -0x1.727569178b6b1p-7, -0x1.862d993d23134p-8,
};

static int check(double actual, double expected)
{
    return actual == expected;
}

static double finite_difference_error(potlibc_handle *handle,
                                      const double geometry[18])
{
    static const double step_factors[] = {
        8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125
    };
    double energy, force[18], displaced[18];
    double worst_best = 0.0;
    if (potlibc_evaluate(handle, 1, geometry, &energy, force) !=
        POTLIBC_SUCCESS) {
        return INFINITY;
    }
    for (size_t coordinate = 0; coordinate < 18; ++coordinate) {
        double best = INFINITY;
        const double base_step =
            cbrt(DBL_EPSILON) * fmax(1.0, fabs(geometry[coordinate]));
        for (size_t step_index = 0;
             step_index < sizeof(step_factors) / sizeof(step_factors[0]);
             ++step_index) {
            double plus, minus;
            const double step = base_step * step_factors[step_index];
            memcpy(displaced, geometry, sizeof(displaced));
            displaced[coordinate] += step;
            if (potlibc_evaluate(handle, 1, displaced, &plus, NULL) !=
                POTLIBC_SUCCESS) {
                return INFINITY;
            }
            displaced[coordinate] -= 2.0 * step;
            if (potlibc_evaluate(handle, 1, displaced, &minus, NULL) !=
                POTLIBC_SUCCESS) {
                return INFINITY;
            }
            const double derivative = (plus - minus) / (2.0 * step);
            const double scale = fmax(1.0, fmax(fabs(derivative),
                                                fabs(force[coordinate])));
            best = fmin(best, fabs(derivative + force[coordinate]) / scale);
        }
        worst_best = fmax(worst_best, best);
    }
    return worst_best;
}

int main(void)
{
    potlibc_handle *handle;
    potlibc_metadata metadata;
    double energy, force[18], product[18], corrected_force[18];
    double corrected_product[18];
    double force_error = 0.0, hvp_error = 0.0;
    double fd_error = 0.0;
    if (potlibc_create("ch52008", &handle) != POTLIBC_SUCCESS) return 1;
    if (potlibc_get_metadata(handle, &metadata) != POTLIBC_SUCCESS) return 2;
    if (metadata.atom_count != 6 ||
        (metadata.capabilities & (POTLIBC_CAP_ENERGY | POTLIBC_CAP_FORCES |
                                  POTLIBC_CAP_HVP)) !=
            (POTLIBC_CAP_ENERGY | POTLIBC_CAP_FORCES | POTLIBC_CAP_HVP))
        return 3;
    if (potlibc_evaluate(handle, 1, x, &energy, force) != POTLIBC_SUCCESS)
        return 4;
    memcpy(corrected_force, force, sizeof(corrected_force));
    if (!check(energy, expected_energy)) {
        fprintf(stderr, "CBE energy got %a expected %a\n", energy,
                expected_energy);
        return 5;
    }
    if (potlibc_hessian_vector_product(handle, 1, x, direction, product) !=
        POTLIBC_SUCCESS) return 7;
    memcpy(corrected_product, product, sizeof(corrected_product));
    static const double geometries[3][18] = {
        {
            2.5448986937601746, 3.0319878740728727, 2.603957609783651,
            1.1100980488026544, 1.3613545004038685, 1.1486133528262976,
            0.013880929393421192, 0.2086334261632169, 2.457516169583444,
            2.4611712350631576, 0.21486419033998002, 0.07550700261319757,
            -0.09665148160819977, 2.4937031206965017, -0.0782634154328256,
            3.3853544918540113, 4.100244654765713, 3.5175358525837166
        },
        {
            2.52, 3.06, 2.58, 1.13, 1.34, 1.17, 0.04, 0.19, 2.49,
            2.43, 0.24, 0.09, -0.07, 2.46, -0.10, 3.41, 4.07, 3.54
        },
        {
            2.57, 3.00, 2.63, 1.08, 1.39, 1.12, -0.02, 0.23, 2.43,
            2.49, 0.18, 0.05, -0.12, 2.52, -0.05, 3.36, 4.13, 3.49
        }
    };
    for (size_t geometry = 0; geometry < 3; ++geometry) {
        fd_error = fmax(fd_error,
                        finite_difference_error(handle, geometries[geometry]));
    }
    const double fd_tolerance =
        32.0 * cbrt(DBL_EPSILON) * cbrt(DBL_EPSILON);
    if (fd_error > fd_tolerance) {
        fprintf(stderr, "CBE corrected gradient FD error %a exceeds %a\n",
                fd_error, fd_tolerance);
        return 6;
    }
    if (potlibc_set_gradient_mode(handle, POTLIBC_GRADIENT_RAW_UPSTREAM) !=
        POTLIBC_SUCCESS) return 8;
    if (potlibc_evaluate(handle, 1, x, &energy, force) != POTLIBC_SUCCESS)
        return 9;
    if (potlibc_hessian_vector_product(handle, 1, x, direction, product) !=
        POTLIBC_SUCCESS) return 10;
    for (size_t i = 0; i < 18; ++i) {
        force_error = fmax(force_error, fabs(force[i] - expected_force[i]));
        if (!check(force[i], expected_force[i])) {
            fprintf(stderr, "CBE raw force[%zu] got %a expected %a\n", i,
                    force[i], expected_force[i]);
            return 11;
        }
        hvp_error = fmax(hvp_error, fabs(product[i] - expected_hvp[i]));
        if (!check(product[i], expected_hvp[i])) {
            fprintf(stderr, "CBE raw HVP[%zu] got %a expected %a\n", i,
                    product[i], expected_hvp[i]);
            return 12;
        }
        if (!check(corrected_force[i], force[i] * CBE_GRADIENT_SCALE) ||
            !check(corrected_product[i], product[i] * CBE_GRADIENT_SCALE)) {
            fprintf(stderr, "CBE correction ratio mismatch at %zu\n", i);
            return 13;
        }
    }
    if (potlibc_set_gradient_mode(
            handle, (potlibc_gradient_mode)-1) !=
        POTLIBC_ERROR_INVALID_GRADIENT_MODE) return 14;
    printf("CBE energy=%a scale=%a fd_best_error=%a raw_force_error=%a "
           "raw_hvp_error=%a\n", energy, CBE_GRADIENT_SCALE, fd_error,
           force_error, hvp_error);
    potlibc_destroy(handle);
    return 0;
}
