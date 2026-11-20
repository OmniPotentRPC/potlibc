/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#define _GNU_SOURCE
#include "potlibc.h"

#include <dlfcn.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef POTLIBC_PLUGIN_DIR
#error "POTLIBC_PLUGIN_DIR must name the installed plugin directory"
#endif

/*
 * Required generated-plugin ABI:
 *   void potlib_plugin_initialize(void);
 *   int potlib_plugin_atom_count(void);  (positive=fixed, zero=variable)
 *   int potlib_plugin_metadata(int atom_count, int *derivatives,
 *       int atomic_numbers[atom_count], double reference[3 * atom_count]);
 *   void potlib_plugin_evaluate(int count, const double *x, double *energy,
 *       double *force);
 * Optional: potlib_plugin_energy_directional, potlib_plugin_hvp, and
 *   void potlib_plugin_gradient_jet(int count, int order, const double *x,
 *       const double *direction, double *jet[order][count][3 * atom_count],
 *       int64_t *operations, double *scale);
 * Metadata returns zero on success and nonzero when atom_count is invalid.
 */
typedef void (*initialize_fn)(void);
typedef int (*atom_count_fn)(void);
typedef int (*metadata_fn)(int, int *, int *, double *);
typedef void (*evaluate_fn)(int, const double *, double *, double *);
typedef void (*directional_fn)(int, const double *, const double *, double *);
typedef void (*hvp_fn)(int, const double *, const double *, double *);
typedef void (*gradient_jet_fn)(int, int, const double *, const double *,
                                double *, int64_t *, double *);
typedef void (*release_fn)(void);

struct potlibc_handle {
    void *image;
    int image_fd;
    char *surface_id;
    int *atomic_numbers;
    double *reference_geometry;
    size_t atom_count;
    uint32_t capabilities;
    evaluate_fn evaluate;
    directional_fn directional;
    hvp_fn hvp;
    gradient_jet_fn gradient_jet;
    release_fn release;
};

static bool valid_id(const char *id)
{
    if (*id == '\0') {
        return false;
    }
    for (const unsigned char *p = (const unsigned char *)id; *p != '\0'; ++p) {
        if (!( (*p >= 'a' && *p <= 'z') || (*p >= '0' && *p <= '9') ||
               *p == '-')) {
            return false;
        }
    }
    return true;
}

static bool finite_values(size_t n, const double *values)
{
    for (size_t i = 0; i < n; ++i) {
        if (!isfinite(values[i])) {
            return false;
        }
    }
    return true;
}

static potlibc_status checked_size(const potlibc_handle *handle, size_t count,
                                   size_t *coordinates)
{
    if (count == 0) {
        return POTLIBC_ERROR_INVALID_COUNT;
    }
    if (count > (size_t)INT_MAX ||
        handle->atom_count > SIZE_MAX / 3 ||
        count > SIZE_MAX / (3 * handle->atom_count) ||
        count * 3 * handle->atom_count > SIZE_MAX / sizeof(double)) {
        return POTLIBC_ERROR_COUNT_RANGE;
    }
    *coordinates = count * 3 * handle->atom_count;
    return POTLIBC_SUCCESS;
}

static void *load_private_image(const char *path, int *retained_fd)
{
    char private_path[64];
    unsigned char buffer[65536];
    int input = -1;
    int image = -1;
    void *handle = NULL;

    input = open(path, O_RDONLY | O_CLOEXEC);
    if (input < 0) {
        return NULL;
    }
    image = memfd_create("potlibc-surface", MFD_CLOEXEC);
    if (image < 0) {
        close(input);
        return NULL;
    }
    for (;;) {
        ssize_t count = read(input, buffer, sizeof(buffer));
        if (count < 0) {
            goto done;
        }
        if (count == 0) {
            break;
        }
        size_t offset = 0;
        while (offset < (size_t)count) {
            ssize_t written =
                write(image, buffer + offset, (size_t)count - offset);
            if (written <= 0) {
                goto done;
            }
            offset += (size_t)written;
        }
    }
    if (snprintf(private_path, sizeof(private_path), "/proc/self/fd/%d",
                 image) < 0) {
        goto done;
    }
    handle = dlopen(private_path, RTLD_NOW | RTLD_LOCAL);
    if (handle != NULL) {
        *retained_fd = image;
        image = -1;
    }
done:
    if (image >= 0) {
        close(image);
    }
    close(input);
    return handle;
}

enum create_mode {
    CREATE_FIXED,
    CREATE_SIZED
};

static potlibc_status create_handle(const char *surface_id,
                                    enum create_mode mode,
                                    size_t requested_atom_count,
                                    potlibc_handle **result)
{
    char path[PATH_MAX];
    const char *directory;
    potlibc_handle *handle;
    initialize_fn initialize;
    atom_count_fn get_atom_count;
    metadata_fn metadata;
    int fixed_atom_count;
    int selected_atom_count;
    int derivatives;
    int length;
    size_t reference_count;

    if (surface_id == NULL || result == NULL) {
        return POTLIBC_ERROR_NULL_POINTER;
    }
    *result = NULL;
    if (!valid_id(surface_id)) {
        return POTLIBC_ERROR_INVALID_SURFACE_ID;
    }
    if (mode == CREATE_SIZED) {
        if (requested_atom_count == 0) {
            return POTLIBC_ERROR_INVALID_ATOM_COUNT;
        }
        if (requested_atom_count > (size_t)INT_MAX) {
            return POTLIBC_ERROR_COUNT_RANGE;
        }
    }
    directory = getenv("POTLIBC_PLUGIN_DIR");
    if (directory == NULL) {
        directory = POTLIBC_PLUGIN_DIR;
    }
    length = snprintf(path, sizeof(path), "%s/libpotlibc-%s.so", directory,
                      surface_id);
    if (length < 0 || (size_t)length >= sizeof(path)) {
        return POTLIBC_ERROR_INVALID_SURFACE_ID;
    }
    handle = calloc(1, sizeof(*handle));
    if (handle == NULL) {
        return POTLIBC_ERROR_ALLOCATION;
    }
    handle->image_fd = -1;
    handle->surface_id = strdup(surface_id);
    if (handle->surface_id == NULL) {
        free(handle);
        return POTLIBC_ERROR_ALLOCATION;
    }
    handle->image = load_private_image(path, &handle->image_fd);
    if (handle->image == NULL) {
        free(handle->surface_id);
        free(handle);
        return POTLIBC_ERROR_SURFACE_UNAVAILABLE;
    }
    initialize = (initialize_fn)dlsym(handle->image, "potlib_plugin_initialize");
    get_atom_count =
        (atom_count_fn)dlsym(handle->image, "potlib_plugin_atom_count");
    metadata = (metadata_fn)dlsym(handle->image, "potlib_plugin_metadata");
    handle->evaluate =
        (evaluate_fn)dlsym(handle->image, "potlib_plugin_evaluate");
    handle->directional =
        (directional_fn)dlsym(handle->image,
                              "potlib_plugin_energy_directional");
    handle->hvp = (hvp_fn)dlsym(handle->image, "potlib_plugin_hvp");
    handle->gradient_jet = (gradient_jet_fn)dlsym(
        handle->image, "potlib_plugin_gradient_jet");
    handle->release =
        (release_fn)dlsym(handle->image, "potlibc_state_release");
    if (initialize == NULL || get_atom_count == NULL || metadata == NULL ||
        handle->evaluate == NULL) {
        potlibc_destroy(handle);
        return POTLIBC_ERROR_PLUGIN_ABI;
    }
    initialize();
    fixed_atom_count = get_atom_count();
    if (fixed_atom_count < 0) {
        potlibc_destroy(handle);
        return POTLIBC_ERROR_PLUGIN_ABI;
    }
    if (fixed_atom_count == 0) {
        if (mode != CREATE_SIZED) {
            potlibc_destroy(handle);
            return POTLIBC_ERROR_ATOM_COUNT_REQUIRED;
        }
        selected_atom_count = (int)requested_atom_count;
    } else {
        if (mode != CREATE_FIXED) {
            potlibc_destroy(handle);
            return POTLIBC_ERROR_ATOM_COUNT_UNEXPECTED;
        }
        selected_atom_count = fixed_atom_count;
    }
    handle->atom_count = (size_t)selected_atom_count;
    if (handle->atom_count > SIZE_MAX / sizeof(*handle->atomic_numbers) ||
        handle->atom_count > SIZE_MAX / 3) {
        potlibc_destroy(handle);
        return POTLIBC_ERROR_COUNT_RANGE;
    }
    reference_count = 3 * handle->atom_count;
    if (reference_count > SIZE_MAX / sizeof(*handle->reference_geometry)) {
        potlibc_destroy(handle);
        return POTLIBC_ERROR_COUNT_RANGE;
    }
    handle->atomic_numbers =
        malloc(handle->atom_count * sizeof(*handle->atomic_numbers));
    handle->reference_geometry =
        malloc(reference_count * sizeof(*handle->reference_geometry));
    if (handle->atomic_numbers == NULL ||
        handle->reference_geometry == NULL) {
        potlibc_destroy(handle);
        return POTLIBC_ERROR_ALLOCATION;
    }
    if (metadata(selected_atom_count, &derivatives, handle->atomic_numbers,
                 handle->reference_geometry) != 0) {
        potlibc_destroy(handle);
        return fixed_atom_count == 0
                   ? POTLIBC_ERROR_INVALID_ATOM_COUNT
                   : POTLIBC_ERROR_PLUGIN_ABI;
    }
    if (derivatives < 0) {
        potlibc_destroy(handle);
        return POTLIBC_ERROR_PLUGIN_ABI;
    }
    handle->capabilities = POTLIBC_CAP_ENERGY;
    if (derivatives > 0 || handle->directional != NULL) {
        handle->capabilities |= POTLIBC_CAP_FORCES;
    }
    if (handle->hvp != NULL) {
        handle->capabilities |= POTLIBC_CAP_HVP;
    }
    if (handle->gradient_jet != NULL) {
        handle->capabilities |= POTLIBC_CAP_GRADIENT_JET;
    }
    *result = handle;
    return POTLIBC_SUCCESS;
}

potlibc_status potlibc_create(const char *surface_id, potlibc_handle **result)
{
    return create_handle(surface_id, CREATE_FIXED, 0, result);
}

potlibc_status potlibc_create_sized(const char *surface_id, size_t atom_count,
                                    potlibc_handle **result)
{
    return create_handle(surface_id, CREATE_SIZED, atom_count, result);
}

void potlibc_destroy(potlibc_handle *handle)
{
    if (handle != NULL) {
        if (handle->image != NULL) {
            if (handle->release != NULL) {
                handle->release();
            }
            dlclose(handle->image);
        }
        if (handle->image_fd >= 0) {
            close(handle->image_fd);
        }
        free(handle->reference_geometry);
        free(handle->atomic_numbers);
        free(handle->surface_id);
        free(handle);
    }
}

potlibc_status potlibc_get_metadata(const potlibc_handle *handle,
                                    potlibc_metadata *metadata)
{
    if (handle == NULL || metadata == NULL) {
        return POTLIBC_ERROR_NULL_POINTER;
    }
    metadata->surface_id = handle->surface_id;
    metadata->atom_count = handle->atom_count;
    metadata->atomic_numbers = handle->atomic_numbers;
    metadata->reference_geometry = handle->reference_geometry;
    metadata->atom_ordering = "POTLIB developer order (atomic numbers array)";
    metadata->length_unit = "bohr";
    metadata->energy_unit = "hartree";
    metadata->force_unit = "hartree/bohr";
    metadata->capabilities = handle->capabilities;
    return POTLIBC_SUCCESS;
}

potlibc_status potlibc_evaluate(potlibc_handle *handle, size_t count,
                                const double *x, double *energy, double *force)
{
    size_t coordinates;
    potlibc_status status;
    double *discarded_force = NULL;

    if (handle == NULL || x == NULL || energy == NULL) {
        return POTLIBC_ERROR_NULL_POINTER;
    }
    if (force != NULL && !(handle->capabilities & POTLIBC_CAP_FORCES)) {
        return POTLIBC_ERROR_UNSUPPORTED;
    }
    status = checked_size(handle, count, &coordinates);
    if (status != POTLIBC_SUCCESS) {
        return status;
    }
    if (!finite_values(coordinates, x)) {
        return POTLIBC_ERROR_NONFINITE_INPUT;
    }
    if (force == NULL) {
        discarded_force = malloc(coordinates * sizeof(*discarded_force));
        if (discarded_force == NULL) {
            return POTLIBC_ERROR_ALLOCATION;
        }
        force = discarded_force;
    }
    handle->evaluate((int)count, x, energy, force);
    if (!finite_values(count, energy) ||
        (!finite_values(coordinates, force) &&
         (handle->capabilities & POTLIBC_CAP_FORCES))) {
        free(discarded_force);
        return POTLIBC_ERROR_NONFINITE_RESULT;
    }
    free(discarded_force);
    return POTLIBC_SUCCESS;
}

potlibc_status potlibc_energy_directional_derivative(
    potlibc_handle *handle, size_t count, const double *x,
    const double *direction, double *derivative)
{
    size_t coordinates;
    potlibc_status status;
    if (handle == NULL || x == NULL || direction == NULL ||
        derivative == NULL) {
        return POTLIBC_ERROR_NULL_POINTER;
    }
    if (handle->directional == NULL) {
        return POTLIBC_ERROR_UNSUPPORTED;
    }
    status = checked_size(handle, count, &coordinates);
    if (status != POTLIBC_SUCCESS) {
        return status;
    }
    if (!finite_values(coordinates, x) ||
        !finite_values(coordinates, direction)) {
        return POTLIBC_ERROR_NONFINITE_INPUT;
    }
    handle->directional((int)count, x, direction, derivative);
    return finite_values(count, derivative)
               ? POTLIBC_SUCCESS
               : POTLIBC_ERROR_NONFINITE_RESULT;
}

potlibc_status potlibc_hessian_vector_product(
    potlibc_handle *handle, size_t count, const double *x,
    const double *direction, double *product)
{
    size_t coordinates;
    potlibc_status status;
    if (handle == NULL || x == NULL || direction == NULL || product == NULL) {
        return POTLIBC_ERROR_NULL_POINTER;
    }
    if (handle->hvp == NULL) {
        return POTLIBC_ERROR_UNSUPPORTED;
    }
    status = checked_size(handle, count, &coordinates);
    if (status != POTLIBC_SUCCESS) {
        return status;
    }
    if (!finite_values(coordinates, x) ||
        !finite_values(coordinates, direction)) {
        return POTLIBC_ERROR_NONFINITE_INPUT;
    }
    handle->hvp((int)count, x, direction, product);
    return finite_values(coordinates, product)
               ? POTLIBC_SUCCESS
               : POTLIBC_ERROR_NONFINITE_RESULT;
}

potlibc_status potlibc_gradient_jet(
    potlibc_handle *handle, size_t count, size_t order, const double *x,
    const double *direction, double *jet, potlibc_jet_report *report)
{
    size_t coordinates;
    potlibc_status status;
    potlibc_jet_report result;
    if (handle == NULL || x == NULL || direction == NULL || jet == NULL) {
        return POTLIBC_ERROR_NULL_POINTER;
    }
    if (handle->gradient_jet == NULL) {
        return POTLIBC_ERROR_UNSUPPORTED;
    }
    if (order < 1 || order > POTLIBC_JET_MAX_ORDER) {
        return POTLIBC_ERROR_JET_ORDER;
    }
    status = checked_size(handle, count, &coordinates);
    if (status != POTLIBC_SUCCESS) {
        return status;
    }
    if (coordinates > SIZE_MAX / (order * sizeof(double))) {
        return POTLIBC_ERROR_COUNT_RANGE;
    }
    if (!finite_values(coordinates, x) ||
        !finite_values(coordinates, direction)) {
        return POTLIBC_ERROR_NONFINITE_INPUT;
    }
    handle->gradient_jet((int)count, (int)order, x, direction, jet,
                         &result.operations, &result.scale);
    if (!finite_values(order * coordinates, jet) || !isfinite(result.scale)) {
        return POTLIBC_ERROR_NONFINITE_RESULT;
    }
    if (report != NULL) {
        *report = result;
    }
    return POTLIBC_SUCCESS;
}

const char *potlibc_status_string(potlibc_status status)
{
    static const char *const messages[] = {
        "success", "null pointer", "invalid surface id", "surface unavailable",
        "plugin ABI mismatch", "allocation failed", "count must be positive",
        "count exceeds implementation range", "non-finite input",
        "non-finite result", "operation unsupported by surface",
        "atom count required for variable-size surface",
        "atom count is only valid for variable-size surfaces",
        "invalid atom count", "jet order outside 1..3"
    };
    size_t index = (size_t)status;
    return index < sizeof(messages) / sizeof(messages[0])
               ? messages[index]
               : "unknown potlibc status";
}
