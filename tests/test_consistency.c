/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include "potlibc.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double dot(size_t n, const double *a, const double *b)
{
    double result = 0.0;
    for (size_t i = 0; i < n; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

static int consistent(size_t operations, double left, double right)
{
    const double rounding =
        sqrt(DBL_EPSILON) * (double)operations * (double)operations;
    const double scale = fabs(left) + fabs(right);
    return fabs(left - right) <= rounding * scale;
}

static int check_surface(const char *test_case)
{
    const char *separator = strrchr(test_case, '=');
    potlibc_handle *handle = NULL;
    potlibc_metadata metadata;
    double energy;
    double derivative;
    double derivative_minus;
    double derivative_plus;
    double derivative_minus_half;
    double derivative_plus_half;
    double finite_difference;
    double finite_difference_half;
    double projected_hvp;
    double roundoff;
    double evaluation_roundoff;
    double step;
    double truncation;
    double *geometry = NULL;
    double *force = NULL;
    double *first = NULL;
    double *second = NULL;
    double *h_first = NULL;
    double *h_second = NULL;
    size_t coordinates;
    size_t selected_atoms;
    char *id = NULL;
    int result = 1;

    if (separator == NULL) {
        fprintf(stderr, "invalid surface test case: %s\n", test_case);
        goto done;
    }
    id = malloc((size_t)(separator - test_case) + 1);
    if (id == NULL) goto done;
    memcpy(id, test_case, (size_t)(separator - test_case));
    id[separator - test_case] = '\0';
    char *end = NULL;
    selected_atoms = strtoull(separator + 1, &end, 10);
    if (end == separator + 1 || *end != '\0') {
        fprintf(stderr, "invalid atom count in test case: %s\n", test_case);
        goto done;
    }
    potlibc_status create_status =
        selected_atoms == 0
            ? potlibc_create(id, &handle)
            : potlibc_create_sized(id, selected_atoms, &handle);
    if (create_status != POTLIBC_SUCCESS ||
        potlibc_get_metadata(handle, &metadata) != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: create or metadata failed\n", id);
        goto done;
    }
    if ((metadata.capabilities &
         (POTLIBC_CAP_ENERGY | POTLIBC_CAP_FORCES | POTLIBC_CAP_HVP)) !=
        (POTLIBC_CAP_ENERGY | POTLIBC_CAP_FORCES | POTLIBC_CAP_HVP)) {
        fprintf(stderr, "%s: required generated capabilities missing\n", id);
        goto done;
    }
    coordinates = metadata.atom_count * 3;
    geometry = malloc(coordinates * sizeof(*geometry));
    force = malloc(coordinates * sizeof(*force));
    first = malloc(coordinates * sizeof(*first));
    second = malloc(coordinates * sizeof(*second));
    h_first = malloc(coordinates * sizeof(*h_first));
    h_second = malloc(coordinates * sizeof(*h_second));
    if (geometry == NULL || force == NULL || first == NULL || second == NULL ||
        h_first == NULL || h_second == NULL) {
        fprintf(stderr, "%s: allocation failed\n", id);
        goto done;
    }
    double reference_scale = 0.0;
    for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        reference_scale = fmax(
            reference_scale, fabs(metadata.reference_geometry[coordinate]));
    }
    for (size_t atom = 0; atom < metadata.atom_count; ++atom) {
        const size_t cell = atom % 8;
        const double spacing = selected_atoms == 0 ? 1.0 : 3.0;
        for (size_t axis = 0; axis < 3; ++axis) {
            const size_t coordinate = 3 * atom + axis;
            if (reference_scale > sqrt(DBL_MIN)) {
                geometry[coordinate] =
                    metadata.reference_geometry[coordinate] +
                    reference_scale * cbrt(DBL_EPSILON) *
                    (double)((coordinate * coordinate) % 13) / 13.0;
            } else {
                const size_t bit = (cell >> axis) & 1;
                const size_t power =
                    axis == 0 ? atom : (axis == 1 ? atom * atom :
                    atom * atom * atom);
                const double divisor =
                    axis == 0 ? 31.0 : (axis == 1 ? 47.0 : 59.0);
                geometry[coordinate] =
                    spacing * (double)(bit + (axis == 0 ? 2 * (atom / 8) : 0)) +
                    (double)(power % (axis == 0 ? 5 : (axis == 1 ? 7 : 11))) /
                    divisor;
            }
        }
    }
    if (strncmp(id, "c2h6", 4) == 0 && metadata.atom_count == 9) {
        static const double ethane_encounter[27] = {
            -1.45, 0.0, 0.0, 1.45, 0.0, 0.0,
            -2.137, 1.944, 0.0, -2.137, -0.972, 1.684,
            -2.137, -0.972, -1.684, 2.137, 1.944, 0.0,
            2.137, -0.972, 1.684, 2.137, -0.972, -1.684,
            0.0, 4.0, 0.0
        };
        memcpy(geometry, ethane_encounter, sizeof(ethane_encounter));
        for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
            geometry[coordinate] +=
                (double)((coordinate * coordinate) % 13) / 997.0;
        }
    }
    if (strncmp(id, "c2h6", 4) == 0 && metadata.atom_count == 10) {
        static const double ethane_oh_encounter[30] = {
            -1.45, 0.0, 0.0, 1.45, 0.0, 0.0,
            -2.137, 1.944, 0.0, -2.137, -0.972, 1.684,
            -2.137, -0.972, -1.684, 2.137, 1.944, 0.0,
            2.137, -0.972, 1.684, 2.137, -0.972, -1.684,
            0.0, 4.0, 0.0, 0.0, 5.8, 0.0
        };
        memcpy(geometry, ethane_oh_encounter, sizeof(ethane_oh_encounter));
        for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
            geometry[coordinate] +=
                (double)((coordinate * coordinate) % 13) / 997.0;
        }
    }
    if (strcmp(id, "phsh-aprp") == 0 && metadata.atom_count == 13) {
        const double ring_radius = 2.64;
        const double hydrogen_radius = ring_radius + 2.06;
        const double pi = acos(-1.0);
        for (size_t atom = 0; atom < 6; ++atom) {
            const double angle = -(double)atom * pi / 3.0;
            geometry[3 * atom] = ring_radius * cos(angle);
            geometry[3 * atom + 1] = ring_radius * sin(angle);
            geometry[3 * atom + 2] = (double)(atom % 2) / 997.0;
        }
        geometry[18] = ring_radius + 3.4;
        geometry[19] = 0.0;
        geometry[20] = 1.0 / 997.0;
        for (size_t atom = 1; atom < 6; ++atom) {
            const double angle = -(double)atom * pi / 3.0;
            const size_t hydrogen = atom + 6;
            geometry[3 * hydrogen] = hydrogen_radius * cos(angle);
            geometry[3 * hydrogen + 1] = hydrogen_radius * sin(angle);
            geometry[3 * hydrogen + 2] =
                (double)((atom * atom) % 7) / 997.0;
        }
        geometry[36] = ring_radius + 5.9;
        geometry[37] = 1.0 / 997.0;
        geometry[38] = -1.0 / 997.0;
    }
    if (strcmp(id, "h3dmbe") == 0 && metadata.atom_count == 3) {
        const double h3_geometry[9] = {
            0.0, 0.0, 0.0, 1.4, 0.0, 0.0,
            0.2, 1.6, 0.1
        };
        memcpy(geometry, h3_geometry, sizeof(h3_geometry));
    }
    if (strncmp(id, "ih2rmc", 6) == 0 && metadata.atom_count == 3) {
        const double ih2_geometry[9] = {
            0.0, 0.0, 0.0, 3.0, 0.1, 0.0,
            4.5, -0.1, 0.05
        };
        memcpy(geometry, ih2_geometry, sizeof(ih2_geometry));
    }
    if ((strcmp(id, "oh2pol") == 0 ||
         strncmp(id, "oh2m3a", 6) == 0) &&
        metadata.atom_count == 3) {
        const double water_geometry[9] = {
            0.0, 0.0, 0.0, 1.8, 0.0, 0.0,
            -0.45, 1.74, 0.03
        };
        memcpy(geometry, water_geometry, sizeof(water_geometry));
    }
    if (strcmp(id, "sioh") == 0 && metadata.atom_count == 3) {
        const double sioh_geometry[9] = {
            0.0, 0.0, 0.0, 3.0, 0.0, 0.0,
            3.9, 1.56, 0.05
        };
        memcpy(geometry, sioh_geometry, sizeof(sioh_geometry));
    }
    for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        first[coordinate] = (double)(coordinate + 1) / (double)coordinates;
        second[coordinate] =
            (double)(coordinates - coordinate) / (double)coordinates;
    }
    potlibc_status evaluate_status =
        potlibc_evaluate(handle, 1, geometry, &energy, force);
    potlibc_status directional_status =
        potlibc_energy_directional_derivative(handle, 1, geometry, first,
                                              &derivative);
    if (evaluate_status != POTLIBC_SUCCESS ||
        directional_status != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: reference evaluation failed (%s, %s)\n", id,
                potlibc_status_string(evaluate_status),
                potlibc_status_string(directional_status));
        goto done;
    }
    if (!consistent(coordinates, derivative, -dot(coordinates, force, first))) {
        fprintf(stderr,
                "%s: energy derivative and gradient disagree (%a, %a)\n",
                id, derivative, -dot(coordinates, force, first));
        goto done;
    }
    if (potlibc_hessian_vector_product(
            handle, 1, geometry, first,
            h_first) != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: HVP evaluation failed\n", id);
        goto done;
    }
    double geometry_scale = 0.0;
    double direction_scale = 0.0;
    for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        geometry_scale = fmax(geometry_scale, fabs(geometry[coordinate]));
        direction_scale = fmax(direction_scale, fabs(second[coordinate]));
    }
    step = cbrt(DBL_EPSILON) *
        fmax(geometry_scale, sqrt(DBL_MIN)) /
        fmax(direction_scale, sqrt(DBL_MIN));
    for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        force[coordinate] =
            geometry[coordinate] + step * second[coordinate];
        h_second[coordinate] =
            geometry[coordinate] - step * second[coordinate];
    }
    if (potlibc_energy_directional_derivative(
            handle, 1, force, first, &derivative_plus) != POTLIBC_SUCCESS ||
        potlibc_energy_directional_derivative(
            handle, 1, h_second, first,
            &derivative_minus) != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: displaced directional evaluation failed\n", id);
        goto done;
    }
    finite_difference =
        (derivative_plus - derivative_minus) / (2.0 * step);
    step *= 0.5;
    for (size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        force[coordinate] =
            geometry[coordinate] + step * second[coordinate];
        h_second[coordinate] =
            geometry[coordinate] - step * second[coordinate];
    }
    if (potlibc_energy_directional_derivative(
            handle, 1, force, first,
            &derivative_plus_half) != POTLIBC_SUCCESS ||
        potlibc_energy_directional_derivative(
            handle, 1, h_second, first,
            &derivative_minus_half) != POTLIBC_SUCCESS) {
        fprintf(stderr, "%s: half-step directional evaluation failed\n", id);
        goto done;
    }
    finite_difference_half =
        (derivative_plus_half - derivative_minus_half) / (2.0 * step);
    truncation =
        fabs((finite_difference_half - finite_difference) / 3.0);
    finite_difference_half +=
        (finite_difference_half - finite_difference) / 3.0;
    projected_hvp = dot(coordinates, second, h_first);
    roundoff = DBL_EPSILON *
        (fabs(derivative_plus_half) + fabs(derivative_minus_half)) / step;
    evaluation_roundoff = sqrt(DBL_EPSILON) *
        (double)coordinates * (double)coordinates * (double)coordinates *
        (fabs(projected_hvp) + fabs(finite_difference_half));
    if (fabs(projected_hvp - finite_difference_half) >
        (double)coordinates *
        (roundoff + truncation) + evaluation_roundoff) {
        fprintf(stderr,
                "%s: HVP disagrees with dual directional derivative "
                "(%a, %a)\n",
                id, projected_hvp, finite_difference_half);
        goto done;
    }
    result = 0;

done:
    free(h_second);
    free(h_first);
    free(second);
    free(first);
    free(force);
    free(geometry);
    potlibc_destroy(handle);
    free(id);
    return result;
}

int main(int argc, char **argv)
{
    for (int surface = 1; surface < argc; ++surface) {
        if (check_surface(argv[surface]) != 0) {
            return surface;
        }
    }
    printf("consistency: %d surfaces passed dual gradient and HVP checks\n",
           argc - 1);
    return 0;
}
