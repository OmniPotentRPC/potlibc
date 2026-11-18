/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Resolves only the symbols required before the gradient-scale ABI and
 * expects a current plugin to refuse that older loader.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>

typedef void (*initialize_fn)(void);
typedef int (*atom_count_fn)(void);
typedef int (*metadata_fn)(int, int *, int *, double *);
typedef void (*evaluate_fn)(int, const double *, double *, double *);

int main(int argc, char **argv)
{
    void *image;
    initialize_fn initialize;
    atom_count_fn atom_count;
    metadata_fn metadata;
    evaluate_fn evaluate;
    int atoms;

    if (argc != 2) {
        fprintf(stderr, "usage: test-plugin-abi-old-loader plugin.so\n");
        return 2;
    }
    image = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (image == NULL) {
        fprintf(stderr, "dlopen: %s\n", dlerror());
        return 3;
    }
    initialize = (initialize_fn)dlsym(image, "potlib_plugin_initialize");
    atom_count = (atom_count_fn)dlsym(image, "potlib_plugin_atom_count");
    metadata = (metadata_fn)dlsym(image, "potlib_plugin_metadata");
    evaluate = (evaluate_fn)dlsym(image, "potlib_plugin_evaluate");
    if (initialize == NULL || atom_count == NULL || metadata == NULL ||
        evaluate == NULL) {
        fprintf(stderr,
                "current plugin is missing a symbol the old loader requires\n");
        dlclose(image);
        return 4;
    }
    initialize();
    atoms = atom_count();
    dlclose(image);
    if (atoms >= 0) {
        fprintf(stderr,
                "old loader accepted a current plugin (atom_count=%d)\n",
                atoms);
        return 5;
    }
    return 0;
}
