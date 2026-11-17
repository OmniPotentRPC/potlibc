/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include "potlibc.h"

#include <stdio.h>
#include <stdlib.h>

static const double geometry[18] = {
    2.5448986937601746, 3.0319878740728727, 2.603957609783651,
    1.1100980488026544, 1.3613545004038685, 1.1486133528262976,
    0.013880929393421192, 0.2086334261632169, 2.457516169583444,
    2.4611712350631576, 0.21486419033998002, 0.07550700261319757,
    -0.09665148160819977, 2.4937031206965017, -0.0782634154328256,
    3.3853544918540113, 4.100244654765713, 3.5175358525837166,
};

int main(void)
{
    potlibc_handle *surface = NULL;
    potlibc_metadata metadata;
    double *direction = NULL;
    double *product = NULL;
    int result = EXIT_FAILURE;

    potlibc_status status = potlibc_create("ch52008", &surface);
    if (status != POTLIBC_SUCCESS) {
        fprintf(stderr, "create: %s\n", potlibc_status_string(status));
        goto done;
    }
    status = potlibc_get_metadata(surface, &metadata);
    if (status != POTLIBC_SUCCESS) {
        fprintf(stderr, "metadata: %s\n", potlibc_status_string(status));
        goto done;
    }
    if ((metadata.capabilities & POTLIBC_CAP_HVP) == 0) {
        fputs("surface does not advertise HVP support\n", stderr);
        goto done;
    }
    if (metadata.atom_count != 6) {
        fputs("unexpected atom count\n", stderr);
        goto done;
    }

    const size_t coordinates = 3 * metadata.atom_count;
    direction = calloc(coordinates, sizeof(*direction));
    product = malloc(coordinates * sizeof(*product));
    if (direction == NULL || product == NULL) {
        fputs("vector allocation failed\n", stderr);
        goto done;
    }
    direction[0] = 1.0;
    status = potlibc_hessian_vector_product(
        surface, 1, geometry, direction, product);
    if (status != POTLIBC_SUCCESS) {
        fprintf(stderr, "HVP: %s\n", potlibc_status_string(status));
        goto done;
    }

    printf("%s: (H e_0)_0 = %.17g hartree/bohr^2\n",
           metadata.surface_id, product[0]);
    result = EXIT_SUCCESS;

done:
    free(product);
    free(direction);
    potlibc_destroy(surface);
    return result;
}
