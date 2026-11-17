/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#ifndef POTLIBC_H
#define POTLIBC_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
# if defined(POTLIBC_BUILDING_LIBRARY)
#  define POTLIBC_API __declspec(dllexport)
# else
#  define POTLIBC_API __declspec(dllimport)
# endif
#elif defined(__GNUC__)
# define POTLIBC_API __attribute__((visibility("default")))
#else
# define POTLIBC_API
#endif

typedef struct potlibc_handle potlibc_handle;

typedef enum potlibc_status {
    POTLIBC_SUCCESS = 0,
    POTLIBC_ERROR_NULL_POINTER,
    POTLIBC_ERROR_INVALID_SURFACE_ID,
    POTLIBC_ERROR_SURFACE_UNAVAILABLE,
    POTLIBC_ERROR_PLUGIN_ABI,
    POTLIBC_ERROR_ALLOCATION,
    POTLIBC_ERROR_INVALID_COUNT,
    POTLIBC_ERROR_COUNT_RANGE,
    POTLIBC_ERROR_NONFINITE_INPUT,
    POTLIBC_ERROR_NONFINITE_RESULT,
    POTLIBC_ERROR_UNSUPPORTED,
    POTLIBC_ERROR_ATOM_COUNT_REQUIRED,
    POTLIBC_ERROR_ATOM_COUNT_UNEXPECTED,
    POTLIBC_ERROR_INVALID_ATOM_COUNT
} potlibc_status;

enum {
    POTLIBC_CAP_ENERGY = UINT32_C(1) << 0,
    POTLIBC_CAP_FORCES = UINT32_C(1) << 1,
    POTLIBC_CAP_HVP = UINT32_C(1) << 2
};

typedef struct potlibc_metadata {
    const char *surface_id;
    size_t atom_count;
    const int *atomic_numbers;
    const double *reference_geometry;
    const char *atom_ordering;
    const char *length_unit;
    const char *energy_unit;
    const char *force_unit;
    uint32_t capabilities;
} potlibc_metadata;

/* Fixed-size surfaces use potlibc_create; variable-size surfaces require
 * potlibc_create_sized. Each function rejects the other surface kind. */
POTLIBC_API potlibc_status potlibc_create(const char *surface_id,
                                          potlibc_handle **handle);
POTLIBC_API potlibc_status potlibc_create_sized(const char *surface_id,
                                                size_t atom_count,
                                                potlibc_handle **handle);
POTLIBC_API void potlibc_destroy(potlibc_handle *handle);
POTLIBC_API potlibc_status potlibc_get_metadata(
    const potlibc_handle *handle, potlibc_metadata *metadata);

/*
 * x and force are row-major [count][3*atom_count]. Coordinates are Bohr,
 * energies Hartree, and forces -dE/dx in Hartree/Bohr. force may be NULL
 * for an energy-only evaluation.
 */
POTLIBC_API potlibc_status potlibc_evaluate(
    potlibc_handle *handle, size_t count, const double *x, double *energy,
    double *force);

/*
 * Forward-mode directional derivative of energy. This exposes the generated
 * dual oracle used to verify each surface's analytic gradient.
 */
POTLIBC_API potlibc_status potlibc_energy_directional_derivative(
    potlibc_handle *handle, size_t count, const double *x,
    const double *direction, double *derivative);

/*
 * Forward-mode analytic H(x)*direction. This is available only when the
 * surface advertises POTLIBC_CAP_HVP.
 */
POTLIBC_API potlibc_status potlibc_hessian_vector_product(
    potlibc_handle *handle, size_t count, const double *x,
    const double *direction, double *product);

POTLIBC_API const char *potlibc_status_string(potlibc_status status);

#ifdef __cplusplus
}
#endif
#endif
