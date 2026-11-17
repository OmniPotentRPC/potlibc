/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#define _POSIX_C_SOURCE 200809L
#include "potlibc.h"

#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef void (*initialize_fn)(void);
typedef int (*atom_count_fn)(void);
typedef int (*metadata_fn)(int, int *, int *, double *);
typedef void (*evaluate_fn)(int, const double *, double *, double *);

enum { sample_count = 9, iterations = 2000 };
static const double cbe_geometry[18] = {
    2.5448986937601746, 3.0319878740728727, 2.603957609783651,
    1.1100980488026544, 1.3613545004038685, 1.1486133528262976,
    0.013880929393421192, 0.2086334261632169, 2.457516169583444,
    2.4611712350631576, 0.21486419033998002, 0.07550700261319757,
    -0.09665148160819977, 2.4937031206965017, -0.0782634154328256,
    3.3853544918540113, 4.100244654765713, 3.5175358525837166,
};

static uint64_t now_ns(void)
{
    struct timespec value;
    if (clock_gettime(CLOCK_MONOTONIC, &value) != 0) {
        return 0;
    }
    return (uint64_t)value.tv_sec * UINT64_C(1000000000) +
           (uint64_t)value.tv_nsec;
}

static void sort(uint64_t *values)
{
    for (size_t i = 1; i < sample_count; ++i) {
        uint64_t value = values[i];
        size_t j = i;
        while (j > 0 && values[j - 1] > value) {
            values[j] = values[j - 1];
            --j;
        }
        values[j] = value;
    }
}

int main(void)
{
    const char *directory = getenv("POTLIBC_PLUGIN_DIR");
    char path[4096];
    void *image;
    initialize_fn initialize;
    atom_count_fn atom_count;
    metadata_fn metadata;
    evaluate_fn evaluate;
    potlibc_handle *handle = NULL;
    uint64_t direct[sample_count], wrapped[sample_count];
    int atoms, derivatives, atomic_numbers[6];
    double geometry[75], energy, force[75];

    if (directory == NULL) return 1;
    if (snprintf(path, sizeof(path), "%s/libpotlibc-ch52008.so", directory) < 0)
        return 2;
    image = dlopen(path, RTLD_NOW | RTLD_LOCAL);
    if (image == NULL) return 3;
    initialize = (initialize_fn)dlsym(image, "potlib_plugin_initialize");
    atom_count = (atom_count_fn)dlsym(image, "potlib_plugin_atom_count");
    metadata = (metadata_fn)dlsym(image, "potlib_plugin_metadata");
    evaluate = (evaluate_fn)dlsym(image, "potlib_plugin_evaluate");
    if (initialize == NULL || atom_count == NULL || metadata == NULL ||
        evaluate == NULL) return 4;
    initialize();
    atoms = atom_count();
    if (atoms != 6) return 5;
    if (metadata(atoms, &derivatives, atomic_numbers, geometry) != 0) return 5;
    for (size_t i = 0; i < 18; ++i) {
        geometry[i] = cbe_geometry[i];
    }
    if (potlibc_create("ch52008", &handle) != POTLIBC_SUCCESS) return 6;

    for (size_t sample = 0; sample < sample_count; ++sample) {
        uint64_t start = now_ns();
        for (int i = 0; i < iterations; ++i) {
            evaluate(1, geometry, &energy, force);
        }
        direct[sample] = now_ns() - start;
        start = now_ns();
        for (int i = 0; i < iterations; ++i) {
            if (potlibc_evaluate(handle, 1, geometry, &energy, force) !=
                POTLIBC_SUCCESS) {
                return 7;
            }
        }
        wrapped[sample] = now_ns() - start;
    }
    sort(direct);
    sort(wrapped);
    potlibc_destroy(handle);
    dlclose(image);
    printf("overhead: n=%d median_direct=%.1f ns median_potlibc=%.1f ns "
           "ratio=%.6f\n",
           sample_count, (double)direct[sample_count / 2] / iterations,
           (double)wrapped[sample_count / 2] / iterations,
           (double)wrapped[sample_count / 2] / direct[sample_count / 2]);
    return 0;
}
