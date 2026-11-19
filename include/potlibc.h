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
    /** The identifier is empty, contains a character other than a
     *  lowercase letter, digit, or hyphen, or does not fit in the plugin
     *  path. A well formed identifier with no loadable plugin returns
     *  POTLIBC_ERROR_SURFACE_UNAVAILABLE. */
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
    POTLIBC_ERROR_INVALID_ATOM_COUNT,
    /** The requested jet order is outside 1..POTLIBC_JET_MAX_ORDER. */
    POTLIBC_ERROR_JET_ORDER
} potlibc_status;

/** Surface capability bits reported by potlibc_metadata.capabilities. */
typedef enum potlibc_capability {
    /** Energy evaluation is available. */
    POTLIBC_CAP_ENERGY = UINT32_C(1) << 0,
    /** Analytic forces are available. */
    POTLIBC_CAP_FORCES = UINT32_C(1) << 1,
    /** Analytic Hessian-vector products are available. */
    POTLIBC_CAP_HVP = UINT32_C(1) << 2,
    /** Directional derivatives of the gradient through third order with a
     *  roundoff report are available (potlibc_gradient_jet). */
    POTLIBC_CAP_GRADIENT_JET = UINT32_C(1) << 3
} potlibc_capability;

/** Highest derivative order of potlibc_gradient_jet. */
#define POTLIBC_JET_MAX_ORDER 3

/**
 * Roundoff report of a gradient jet.
 *
 * @c operations counts the forward-mode operations of the whole batch and
 * @c scale is the largest magnitude among the gradient and the returned
 * derivative orders. Their product with DBL_EPSILON is a first-order estimate
 * of the roundoff of the derivatives from an operation count; it is not a
 * rigorous bound.
 */
typedef struct potlibc_jet_report {
    int64_t operations;
    double scale;
} potlibc_jet_report;

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
 * A zero atom count returns POTLIBC_ERROR_INVALID_ATOM_COUNT.
 * Fixed-size surfaces return POTLIBC_ERROR_ATOM_COUNT_UNEXPECTED.
 * A variable-size surface whose metadata rejects the count returns
 * POTLIBC_ERROR_INVALID_ATOM_COUNT.
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
 * Evaluate directional derivatives of the analytic gradient through a given
 * order.
 *
 * @param handle Live surface handle.
 * @param count Number of geometries.
 * @param order Highest derivative order, in 1..POTLIBC_JET_MAX_ORDER.
 * @param x Row-major coordinates with shape [count][3 * atom_count], in Bohr.
 * @param direction Row-major direction vectors with the same shape as @p x.
 * @param jet Row-major output with shape [order][count][3 * atom_count].
 * @param report Optional roundoff report, written only on success.
 * @return POTLIBC_SUCCESS, POTLIBC_ERROR_UNSUPPORTED,
 * POTLIBC_ERROR_JET_ORDER, or another diagnostic status.
 *
 * jet[k - 1] is the k-th derivative with respect to t of the gradient of
 * E(x + t direction) at t = 0, in Hartree/Bohr^(k + 1) times the direction
 * to the power k. Order 1 is the Hessian-vector product, order 2 contracts the
 * third derivative tensor with the direction twice and order 3 contracts the
 * fourth derivative tensor with it three times. They are exact forward-mode
 * values, not Taylor coefficients, and no differencing step is used.
 *
 * The operation is available when POTLIBC_CAP_GRADIENT_JET is set.
 */
POTLIBC_API potlibc_status potlibc_gradient_jet(
    potlibc_handle *handle, size_t count, size_t order, const double *x,
    const double *direction, double *jet, potlibc_jet_report *report);

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
