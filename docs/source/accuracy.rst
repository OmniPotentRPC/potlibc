.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

Accuracy and thread safety
==========================

What the HVP computes
---------------------

Generated mixed dual values carry a coordinate seed, a direction seed, and
their mixed term. Forces read the coordinate derivative; HVPs read the mixed
derivative. The implementation therefore computes the derivative of the
analytic gradient in one pass and never forms a finite-difference Hessian.

Branches are selected from primal values. At nonsmooth points, ``abs(0)`` and
``sqrt(0)`` use a zero subgradient, clipping ties preserve a constant bound,
and integral real powers use integer arithmetic for negative bases. An HVP at
a branch or clipping boundary is the derivative of that selected convention,
not a claim that a classical Hessian exists there.

Measured checks
---------------

The release test suite reports:

.. code-block:: text

   CBE energy=0x1.885b2451b0c6p-6 energy_abs_error=0x0p+0 force_max_abs_error=0x0p+0 hvp_max_abs_error=0x0p+0
   threading: 128 concurrent handles across 64 threads and 6 surfaces, bit-identical energy and forces

The CBE regression compares one energy, all 18 force components, and all 18
HVP components bit-for-bit with committed hexadecimal references. Its maximum
absolute error is exactly zero for each quantity.

The complete Meson run passed 400 of 400 tests. This includes 387 of 387
per-surface dual consistency cases. For every catalog surface, that test:

#. compares the dual directional energy derivative with
   ``-dot(force, direction)``;
#. compares a projected analytic HVP with a Richardson-extrapolated centered
   difference of dual directional derivatives; and
#. scales the acceptance bound from double-precision roundoff, measured
   truncation, coordinate dimension, and displaced derivative magnitudes.

The finite-difference probe is an independent check, not the production HVP
algorithm. The adaptive bound is preferable to one global decimal tolerance
because the catalog spans dimensions and energy scales that differ by orders
of magnitude.

Concurrency contract
--------------------

Legacy surfaces retain mutable Fortran ``COMMON`` and ``SAVE`` state.
``potlibc`` gives each handle a separately loaded plugin image. Adapters
other than CBE (``ch52008``) restore that image's initialized writable-data
snapshot before each seed evaluation. CBE does not restore a snapshot.
No global evaluation lock is used.

The parallel regression holds 128 handles live across 64 threads for six
surface families and requires energy and every force bit to match serial
evaluation. The same case is available under ThreadSanitizer through the
``test-tsan`` target when the compiler accepts ``-fsanitize=thread`` at
compile time and at link time.

The supported contract is concurrent calls on distinct handles. Calls on one
handle must be externally serialized.
