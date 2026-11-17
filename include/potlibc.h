/* SPDX-FileCopyrightText: 2026 Rohit Goswami
 * SPDX-License-Identifier: BSD-3-Clause
 */
#ifndef POTLIBC_H
#define POTLIBC_H

#include <stddef.h>
#include <stdint.h>

/**
 * @file potlibc.h
 * @brief Stable C ABI for the PotLib 2023 surface catalog.
 *
 * Coordinates are in Bohr, energies are in Hartree, forces are in
 * Hartree/Bohr, and Hessian-vector products are in Hartree/Bohr^2.
 * Arrays are contiguous and row-major. Unless a function says otherwise,
 * pointers must remain valid only for the duration of the call.
 */

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

/** Opaque, independently isolated surface instance. */
typedef struct potlibc_handle potlibc_handle;

/** Status returned by every fallible operation. */
typedef enum potlibc_status {
    /** The operation completed successfully. */
    POTLIBC_SUCCESS = 0,
    /** A required pointer argument was NULL. */
    POTLIBC_ERROR_NULL_POINTER,
    /** No catalog entry has the requested surface identifier. */
    POTLIBC_ERROR_INVALID_SURFACE_ID,
    /** The surface plugin could not be opened. */
    POTLIBC_ERROR_SURFACE_UNAVAILABLE,
    /** The loaded plugin does not implement the expected internal ABI. */
    POTLIBC_ERROR_PLUGIN_ABI,
    /** Memory allocation failed. */
    POTLIBC_ERROR_ALLOCATION,
    /** A batch count is zero. */
    POTLIBC_ERROR_INVALID_COUNT,
    /** A count or buffer size exceeds the supported integer range. */
    POTLIBC_ERROR_COUNT_RANGE,
    /** At least one coordinate or direction component is not finite. */
    POTLIBC_ERROR_NONFINITE_INPUT,
    /** The surface returned a non-finite energy, force, or product. */
    POTLIBC_ERROR_NONFINITE_RESULT,
    /** The requested operation is not advertised by this surface. */
    POTLIBC_ERROR_UNSUPPORTED,
    /** A variable-size surface was created without an atom count. */
    POTLIBC_ERROR_ATOM_COUNT_REQUIRED,
    /** A fixed-size surface was created with an atom count. */
    POTLIBC_ERROR_ATOM_COUNT_UNEXPECTED,
    /** The requested atom count is not valid for the surface. */
    POTLIBC_ERROR_INVALID_ATOM_COUNT
} potlibc_status;

/** Surface capability bits reported by potlibc_metadata.capabilities. */
typedef enum potlibc_capability {
    /** Energy evaluation is available. */
    POTLIBC_CAP_ENERGY = UINT32_C(1) << 0,
    /** Analytic forces are available. */
    POTLIBC_CAP_FORCES = UINT32_C(1) << 1,
    /** Analytic Hessian-vector products are available. */
    POTLIBC_CAP_HVP = UINT32_C(1) << 2
} potlibc_capability;

/**
 * Immutable metadata owned by a potlibc_handle.
 *
 * String and array pointers remain valid until the associated handle is
 * destroyed. Callers must not free or modify them.
 */
typedef struct potlibc_metadata {
    /** Canonical catalog identifier accepted by the creation functions. */
    const char *surface_id;
    /** Number of atoms expected in each geometry. */
    size_t atom_count;
    /** Atomic numbers in the deposited atom ordering; length atom_count. */
    const int *atomic_numbers;
    /** Cartesian reference geometry; length 3 * atom_count. */
    const double *reference_geometry;
    /** Human-readable description of the deposited atom ordering. */
    const char *atom_ordering;
    /** Coordinate unit, currently "bohr". */
    const char *length_unit;
    /** Energy unit, currently "hartree". */
    const char *energy_unit;
    /** Force unit, currently "hartree/bohr". */
    const char *force_unit;
    /** Bitwise OR of POTLIBC_CAP_* values. */
    uint32_t capabilities;
} potlibc_metadata;

/**
 * Create an isolated handle for a fixed-size surface.
 *
 * @param surface_id NUL-terminated canonical catalog identifier.
 * @param handle Receives the new handle on success and NULL on failure.
 * @return POTLIBC_SUCCESS or a diagnostic status.
 *
 * Variable-size surfaces return POTLIBC_ERROR_ATOM_COUNT_REQUIRED.
 * The returned handle may be used concurrently with other handles, but calls
 * on the same handle must not overlap.
 */
POTLIBC_API potlibc_status potlibc_create(const char *surface_id,
                                          potlibc_handle **handle);

/**
 * Create an isolated handle for a variable-size surface.
 *
 * @param surface_id NUL-terminated canonical catalog identifier.
 * @param atom_count Number of atoms in every geometry evaluated by the handle.
 * @param handle Receives the new handle on success and NULL on failure.
 * @return POTLIBC_SUCCESS or a diagnostic status.
 *
 * Fixed-size surfaces return POTLIBC_ERROR_ATOM_COUNT_UNEXPECTED.
 */
POTLIBC_API potlibc_status potlibc_create_sized(const char *surface_id,
                                                size_t atom_count,
                                                potlibc_handle **handle);

/**
 * Destroy a handle. Passing NULL is permitted.
 *
 * The caller must ensure no operation on @p handle is in progress.
 */
POTLIBC_API void potlibc_destroy(potlibc_handle *handle);

/**
 * Read immutable metadata for a handle.
 *
 * @param handle Live surface handle.
 * @param metadata Receives borrowed metadata pointers.
 * @return POTLIBC_SUCCESS or POTLIBC_ERROR_NULL_POINTER.
 */
POTLIBC_API potlibc_status potlibc_get_metadata(
    const potlibc_handle *handle, potlibc_metadata *metadata);

/**
 * Evaluate a batch of geometries.
 *
 * @param handle Live surface handle.
 * @param count Number of geometries; must be greater than zero.
 * @param x Row-major coordinates with shape [count][3 * atom_count], in Bohr.
 * @param energy Output array of length count, in Hartree.
 * @param force Optional row-major output with the same shape as @p x,
 * containing -dE/dx in Hartree/Bohr.
 * @return POTLIBC_SUCCESS or a diagnostic status.
 */
POTLIBC_API potlibc_status potlibc_evaluate(
    potlibc_handle *handle, size_t count, const double *x, double *energy,
    double *force);

/**
 * Evaluate the forward-mode directional derivative of the energy.
 *
 * @param handle Live surface handle.
 * @param count Number of geometries.
 * @param x Row-major coordinates with shape [count][3 * atom_count].
 * @param direction Row-major direction vectors with the same shape as @p x.
 * @param derivative Output array of length count containing grad(E) dot
 * direction in Hartree/Bohr.
 * @return POTLIBC_SUCCESS or a diagnostic status.
 *
 * This independent dual-number path is also the oracle used to verify each
 * generated analytic force adapter.
 */
POTLIBC_API potlibc_status potlibc_energy_directional_derivative(
    potlibc_handle *handle, size_t count, const double *x,
    const double *direction, double *derivative);

/**
 * Evaluate an analytic Hessian-vector product without forming a Hessian.
 *
 * @param handle Live surface handle.
 * @param count Number of geometries.
 * @param x Row-major coordinates with shape [count][3 * atom_count].
 * @param direction Row-major direction vectors with the same shape as @p x.
 * @param product Row-major output with the same shape as @p x, in
 * Hartree/Bohr^2.
 * @return POTLIBC_SUCCESS, POTLIBC_ERROR_UNSUPPORTED, or another diagnostic
 * status.
 *
 * The operation is available when POTLIBC_CAP_HVP is set. It propagates a
 * mixed forward-dual seed through the generated surface and never allocates a
 * full Hessian.
 */
POTLIBC_API potlibc_status potlibc_hessian_vector_product(
    potlibc_handle *handle, size_t count, const double *x,
    const double *direction, double *product);

/**
 * Return a static English description of a status value.
 *
 * The returned pointer is never NULL and must not be freed.
 */
POTLIBC_API const char *potlibc_status_string(potlibc_status status);

#ifdef __cplusplus
}
#endif
#endif
