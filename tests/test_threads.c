/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#define _GNU_SOURCE
#include "potlibc.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { surface_count = 6, thread_count = 64, handles_per_thread = 2 };
static const char *const surfaces[surface_count] = {
    "ch52008", "brh2sec2001", "fh2sk2001",
    "h3tk2001", "ih2raff2001", "clhclbcmr2001",
};
static const double cbe_geometry[18] = {
    2.5448986937601746, 3.0319878740728727, 2.603957609783651,
    1.1100980488026544, 1.3613545004038685, 1.1486133528262976,
    0.013880929393421192, 0.2086334261632169, 2.457516169583444,
    2.4611712350631576, 0.21486419033998002, 0.07550700261319757,
    -0.09665148160819977, 2.4937031206965017, -0.0782634154328256,
    3.3853544918540113, 4.100244654765713, 3.5175358525837166,
};
static const double triatomic_geometry[9] = {
    0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 1.7, 0.0,
};
struct reference {
    size_t coordinates;
    double geometry[75];
    double energy;
    double force[75];
};
struct worker_data {
    size_t index;
};

static pthread_barrier_t barrier;
static atomic_int failed;
static struct reference serial[surface_count];

static void *worker(void *argument)
{
    const struct worker_data *data = argument;
    potlibc_handle *handles[handles_per_thread] = {NULL, NULL};
    size_t selected[handles_per_thread];
    double energy[handles_per_thread];
    double force[handles_per_thread][75];

    for (size_t i = 0; i < handles_per_thread; ++i) {
        selected[i] = (data->index * handles_per_thread + i) % surface_count;
        if (potlibc_create(surfaces[selected[i]], &handles[i]) !=
            POTLIBC_SUCCESS) {
            atomic_store(&failed, 1);
        }
    }
    pthread_barrier_wait(&barrier);
    for (size_t i = 0; i < handles_per_thread; ++i) {
        const struct reference *expected = &serial[selected[i]];
        if (handles[i] != NULL &&
            (potlibc_evaluate(handles[i], 1, expected->geometry, &energy[i],
                              force[i]) !=
                POTLIBC_SUCCESS ||
             memcmp(&energy[i], &expected->energy, sizeof(energy[i])) != 0 ||
             memcmp(force[i], expected->force,
                    expected->coordinates * sizeof(force[i][0])) != 0)) {
            atomic_store(&failed, 1);
        }
    }
    pthread_barrier_wait(&barrier);
    for (size_t i = 0; i < handles_per_thread; ++i) {
        potlibc_destroy(handles[i]);
    }
    return NULL;
}

int main(void)
{
    pthread_t threads[thread_count];
    struct worker_data data[thread_count];

    for (size_t i = 0; i < surface_count; ++i) {
        potlibc_handle *handle = NULL;
        potlibc_metadata metadata;
        if (potlibc_create(surfaces[i], &handle) != POTLIBC_SUCCESS) return 1;
        if (potlibc_get_metadata(handle, &metadata) != POTLIBC_SUCCESS) return 2;
        serial[i].coordinates = 3 * metadata.atom_count;
        if (i == 0 && serial[i].coordinates == 18) {
            memcpy(serial[i].geometry, cbe_geometry, sizeof(cbe_geometry));
        } else if (i != 0 && serial[i].coordinates == 9) {
            memcpy(serial[i].geometry, triatomic_geometry,
                   sizeof(triatomic_geometry));
        } else {
            return 3;
        }
        if (potlibc_evaluate(handle, 1, serial[i].geometry, &serial[i].energy,
                             serial[i].force) != POTLIBC_SUCCESS) return 4;
        potlibc_destroy(handle);
    }
    if (pthread_barrier_init(&barrier, NULL, thread_count) != 0) return 4;
    for (size_t i = 0; i < thread_count; ++i) {
        data[i].index = i;
        if (pthread_create(&threads[i], NULL, worker, &data[i]) != 0) return 5;
    }
    for (size_t i = 0; i < thread_count; ++i) {
        if (pthread_join(threads[i], NULL) != 0) return 6;
    }
    pthread_barrier_destroy(&barrier);
    printf("threading: %d concurrent handles across %d threads and %d surfaces, "
           "bit-identical energy and forces\n",
           thread_count * handles_per_thread, thread_count, surface_count);
    return atomic_load(&failed);
}
