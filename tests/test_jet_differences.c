/* SPDX-FileCopyrightText: 2026 Rohit Goswami */
/* SPDX-License-Identifier: BSD-3-Clause */

/* Independent oracle for the gradient jet.
 *
 * The jet comes from forward-mode arithmetic on the transformed potential
 * source. This test does not use that arithmetic: it differentiates the
 * lower-order quantities along the same line by central differences with one
 * Richardson step and compares.
 *
 *   g(t)  = gradient of E at x + t v, from potlibc_evaluate (the force, negated)
 *   j1(t) = order-1 jet at x + t v with direction v, which equals g'(t)
 *
 * so the jet at x must satisfy
 *   order 1 = g'(0),   order 2 = j1'(0) = g''(0),   order 3 = j1''(0) = g'''(0).
 *
 * Central differences with step h have error c h^2; combining the steps h and
 * h/2 as (4 D(h/2) - D(h)) / 3 removes that term and leaves c' h^4. The
 * truncation error of the combined value is estimated as the distance between
 * the combined values at h and at 2h (about 15 times the error at h). The
 * roundoff error of a k-th difference is eps * M / h^k with M the largest
 * magnitude among the sampled values, times 16 for the sums involved. A jet
 * component passes when it is within the truncation estimate plus the roundoff
 * estimate of the difference value, and the test also fails when that
 * tolerance is not below 1e-6 of the largest jet component, so a loose oracle
 * cannot pass. The geometries are the saddle and a displaced point of the jet
 * goldens; the differences sample only x + t v with |t| <= 2 h, and agreement
 * of the combined values across h confirms that no branch of the potential
 * switches inside that range.
 *
 * Negative control: the same check on a copy of the jet whose first order-3
 * component is raised by a relative 1e-5 must fail at both geometries, so the
 * oracle is shown to resolve an error far below the size of the jet. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "potlibc.h"
#include "ch52008_jet_oracle.h"

#define N 18

static potlibc_handle *handle;

static int quiet;

static int line_gradient(const double *x, const double *v, double t, double *g,
                         double *j1)
{
    double point[N];
    double energy;
    double force[N];
    double jet[N];
    size_t i;
    for (i = 0; i < N; i++) {
        point[i] = x[i] + t * v[i];
    }
    if (potlibc_evaluate(handle, 1, point, &energy, force) != POTLIBC_SUCCESS ||
        potlibc_gradient_jet(handle, 1, 1, point, v, jet, NULL) !=
            POTLIBC_SUCCESS) {
        return 1;
    }
    for (i = 0; i < N; i++) {
        g[i] = -force[i];
        j1[i] = jet[i];
    }
    return 0;
}

typedef struct {
    double g[N];
    double j1[N];
} sample;

static int take(const double *x, const double *v, double t, sample *s)
{
    return line_gradient(x, v, t, s->g, s->j1);
}

/* Combined (Richardson) k-th derivative of f(t) at 0 from samples at 0, +-h,
 * +-h/2, where f is g (kind 0) or j1 (kind 1); first(order 1) or second. */
static double combined(const sample *at0, const sample *p, const sample *m,
                       const sample *ph, const sample *mh, double h, int kind,
                       int derivative, size_t i, double *roundoff)
{
    const double f0 = kind ? at0->j1[i] : at0->g[i];
    const double fp = kind ? p->j1[i] : p->g[i];
    const double fm = kind ? m->j1[i] : m->g[i];
    const double fph = kind ? ph->j1[i] : ph->g[i];
    const double fmh = kind ? mh->j1[i] : mh->g[i];
    double big = fmax(fmax(fabs(f0), fabs(fp)), fmax(fabs(fm), fmax(fabs(fph), fabs(fmh))));
    double d_h;
    double d_half;
    if (derivative == 1) {
        d_h = (fp - fm) / (2.0 * h);
        d_half = (fph - fmh) / h;
        *roundoff = 16.0 * DBL_EPSILON * big / (h / 2.0);
    } else {
        d_h = (fp - 2.0 * f0 + fm) / (h * h);
        d_half = (fph - 2.0 * f0 + fmh) / ((h / 2.0) * (h / 2.0));
        *roundoff = 16.0 * DBL_EPSILON * big / ((h / 2.0) * (h / 2.0));
    }
    return (4.0 * d_half - d_h) / 3.0;
}

static int check_order(const char *label, const double *x, const double *v,
                       const double *jet_block, int kind, int derivative,
                       double h)
{
    sample at0, p, m, ph, mh, p2, m2;
    double largest = 0.0;
    int failure = 0;
    size_t i;
    if (take(x, v, 0.0, &at0) || take(x, v, h, &p) || take(x, v, -h, &m) ||
        take(x, v, h / 2.0, &ph) || take(x, v, -h / 2.0, &mh) ||
        take(x, v, 2.0 * h, &p2) || take(x, v, -2.0 * h, &m2)) {
        return 1;
    }
    for (i = 0; i < N; i++) {
        largest = fmax(largest, fabs(jet_block[i]));
    }
    for (i = 0; i < N; i++) {
        double round_h;
        double round_2h;
        double at_h = combined(&at0, &p, &m, &ph, &mh, h, kind, derivative, i, &round_h);
        double at_2h = combined(&at0, &p2, &m2, &p, &m, 2.0 * h, kind, derivative, i, &round_2h);
        double tolerance = fabs(at_h - at_2h) + round_h;
        double gap = fabs(jet_block[i] - at_h);
        if (tolerance > 1e-6 * largest) {
            if (!quiet) fprintf(stderr, "%s[%zu]: oracle tolerance %.3e too loose against %.3e\n",
                    label, i, tolerance, largest);
            failure = 1;
        }
        if (gap > tolerance) {
            if (!quiet) fprintf(stderr, "%s[%zu]: jet %.17g, differences %.17g, gap %.3e > %.3e\n",
                    label, i, jet_block[i], at_h, gap, tolerance);
            failure = 1;
        }
    }
    return failure;
}

static int check_point(const char *name, const double *x, const double *v,
                       double order3_bump)
{
    double jet[3 * N];
    double h = 0.02;
    int failure = 0;
    char label[64];
    if (potlibc_gradient_jet(handle, 1, 3, x, v, jet, NULL) != POTLIBC_SUCCESS) {
        return 1;
    }
    jet[2 * N] *= 1.0 + order3_bump;
    snprintf(label, sizeof label, "%s order 1 vs g'", name);
    failure |= check_order(label, x, v, jet, 0, 1, h);
    snprintf(label, sizeof label, "%s order 2 vs j1'", name);
    failure |= check_order(label, x, v, jet + N, 1, 1, h);
    snprintf(label, sizeof label, "%s order 3 vs j1''", name);
    failure |= check_order(label, x, v, jet + 2 * N, 1, 2, h);
    return failure;
}

int main(void)
{
    int failure;
    if (potlibc_create("ch52008", &handle) != POTLIBC_SUCCESS) {
        return 2;
    }
    failure = check_point("saddle", saddle_x, saddle_v, 0.0);
    failure |= check_point("displaced", displaced_x, displaced_v, 0.0);
    quiet = 1;
    if (check_point("saddle bumped", saddle_x, saddle_v, 1e-5) == 0 ||
        check_point("displaced bumped", displaced_x, displaced_v, 1e-5) == 0) {
        fprintf(stderr, "a 1e-5 relative bump of an order-3 component was not detected\n");
        failure = 1;
    }
    potlibc_destroy(handle);
    return failure;
}
