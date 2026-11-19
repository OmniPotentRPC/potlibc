/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include "potlibc.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

enum { coordinates = 9 };

static int evaluate(
    potlibc_handle *handle, const double x[coordinates],
    double force[coordinates])
{
    double energy;

    return potlibc_evaluate(handle, 1, x, &energy, force) == POTLIBC_SUCCESS;
}

static void collinear(double separation, double x[coordinates])
{
    memset(x, 0, coordinates * sizeof(*x));
    x[0] = -separation;
    x[6] = separation;
}

static int find_lsth_saddle(potlibc_handle *handle, double *separation)
{
    const double step = 1.0e-5;
    double x[coordinates];
    double plus[coordinates];
    double minus[coordinates];
    double force[coordinates];
    double force_plus[coordinates];
    double force_minus[coordinates];

    *separation = 1.75;
    for (int iteration = 0; iteration < 12; ++iteration) {
        collinear(*separation, x);
        collinear(*separation + step, plus);
        collinear(*separation - step, minus);
        if (!evaluate(handle, x, force) ||
            !evaluate(handle, plus, force_plus) ||
            !evaluate(handle, minus, force_minus)) {
            return 0;
        }
        const double derivative = -force[0] + force[6];
        const double curvature =
            ((-force_plus[0] + force_plus[6]) -
             (-force_minus[0] + force_minus[6])) / (2.0 * step);
        if (fabs(derivative) < 1.0e-11) return 1;
        if (curvature == 0.0 || !isfinite(curvature)) return 0;
        *separation -= derivative / curvature;
    }
    return 0;
}

static int check_hessian(
    const char *surface, potlibc_handle *handle, const double x[coordinates])
{
    const double step = 1.0e-6;
    double direction[coordinates] = {0.0};
    double product[coordinates];
    double plus[coordinates];
    double minus[coordinates];
    double force_plus[coordinates];
    double force_minus[coordinates];
    double max_error = 0.0;
    double max_reference = 0.0;

    for (int column = 0; column < coordinates; ++column) {
        memset(direction, 0, sizeof(direction));
        direction[column] = 1.0;
        if (potlibc_hessian_vector_product(
                handle, 1, x, direction, product) != POTLIBC_SUCCESS) {
            fprintf(stderr, "%s: HVP evaluation failed\n", surface);
            return 0;
        }
        memcpy(plus, x, sizeof(plus));
        memcpy(minus, x, sizeof(minus));
        plus[column] += step;
        minus[column] -= step;
        if (!evaluate(handle, plus, force_plus) ||
            !evaluate(handle, minus, force_minus)) {
            fprintf(stderr, "%s: displaced gradient evaluation failed\n",
                    surface);
            return 0;
        }
        for (int row = 0; row < coordinates; ++row) {
            const double finite_difference =
                -(force_plus[row] - force_minus[row]) / (2.0 * step);
            const double error = fabs(product[row] - finite_difference);
            max_error = fmax(max_error, error);
            max_reference = fmax(max_reference, fabs(finite_difference));
        }
    }
    printf("%s: maximum Hessian error %.16e (reference %.16e)\n",
           surface, max_error, max_reference);
    if (max_error > 1.0e-7 * fmax(1.0, max_reference)) {
        fprintf(stderr, "%s: dual2 Hessian disagrees with gradient FD\n",
                surface);
        return 0;
    }
    return 1;
}

static int check_surface(
    const char *surface, const char *label, const double x[coordinates])
{
    potlibc_handle *handle = NULL;
    int result = 0;
    potlibc_status status;

    status = potlibc_create(surface, &handle);
    if (status != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: handle creation failed: %s\n", surface,
                potlibc_status_string(status));
        return 0;
    }
    result = check_hessian(label, handle, x);
    potlibc_destroy(handle);
    return result;
}

static int print_generic_value(
    const char *surface, const double x[coordinates])
{
    potlibc_handle *handle = NULL;
    double force[coordinates];
    double energy;

    if (potlibc_create(surface, &handle) != POTLIBC_SUCCESS ||
        potlibc_evaluate(handle, 1, x, &energy, force) != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: generic evaluation failed\n", surface);
        potlibc_destroy(handle);
        return 0;
    }
    printf("%s generic energy %a forces", surface, energy);
    for (int coordinate = 0; coordinate < coordinates; ++coordinate) {
        printf(" %a", force[coordinate]);
    }
    putchar('\n');
    potlibc_destroy(handle);
    return 1;
}

int main(int argc, char **argv)
{
    potlibc_handle *lsth = NULL;
    double saddle[coordinates];
    double symmetric[coordinates];
    const double generic[coordinates] = {
        -1.0, 0.2, 0.1, 0.3, -0.4, 1.7, 2.1, 0.5, -0.2
    };
    double separation;

    if (argc == 2 && strcmp(argv[1], "chc") == 0) {
        collinear(1.70, symmetric);
        printf("chc symmetric separation %.16e\n", 1.70);
        return check_surface("chc", "chc symmetric", symmetric) ? 0 : 1;
    }
    if (argc != 1) {
        fprintf(stderr, "unknown symmetric Hessian surface\n");
        return 1;
    }
    if (!print_generic_value("h3lsth", generic)) {
        return 1;
    }
    if (potlibc_create("h3lsth", &lsth) != POTLIBC_SUCCESS ||
        !find_lsth_saddle(lsth, &separation)) {
        fprintf(stderr, "h3lsth: failed to locate collinear saddle\n");
        potlibc_destroy(lsth);
        return 1;
    }
    collinear(separation, saddle);
    printf("h3lsth saddle separation %.16e\n", separation);
    if (!check_hessian("h3lsth saddle", lsth, saddle)) {
        potlibc_destroy(lsth);
        return 1;
    }
    memcpy(symmetric, saddle, sizeof(symmetric));
    symmetric[1] = 0.37;
    symmetric[4] = 0.37;
    symmetric[7] = 0.37;
    if (!check_hessian("h3lsth symmetric", lsth, symmetric)) {
        potlibc_destroy(lsth);
        return 1;
    }
    potlibc_destroy(lsth);
    return 0;
}
