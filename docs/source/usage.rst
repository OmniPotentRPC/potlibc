.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

C usage
=======

Create, inspect, evaluate, and destroy
--------------------------------------

The catalog identifier selects a surface. Fixed-size surfaces use
:c:func:`potlibc_create`; variable-size surfaces require
:c:func:`potlibc_create_sized`. Metadata pointers are borrowed from the handle
and remain valid until :c:func:`potlibc_destroy`.

This complete example evaluates one six-atom geometry:

.. literalinclude:: ../../examples/energy.c
   :language: c
   :linenos:

``force`` is :math:`-\nabla E`. Pass ``NULL`` for it when only energy is
required. For a batch of ``count`` geometries, flatten coordinates and forces
as ``[count][atom_count][3]`` and provide ``count`` energy values.

Hessian-vector products
-----------------------

:c:func:`potlibc_hessian_vector_product` computes
:math:`H(x)v = \nabla^2 E(x)v` directly with mixed forward dual numbers. It
does not allocate or expose a dense Hessian.

.. literalinclude:: ../../examples/hvp.c
   :language: c
   :linenos:

Gradient jets
-------------

:c:func:`potlibc_gradient_jet` evaluates directional derivatives of the
analytic gradient through a chosen order. For a unit direction :math:`v`, order
``k`` is the :math:`k`-th derivative with respect to :math:`t` of
:math:`\nabla E(x + t v)` at :math:`t = 0`. Order 1 matches
:c:func:`potlibc_hessian_vector_product`; higher orders contract the energy
tensor along :math:`v` without finite differencing.

Only surfaces that advertise :c:macro:`POTLIBC_CAP_GRADIENT_JET` support the
call. ``order`` must lie between 1 and :c:macro:`POTLIBC_JET_MAX_ORDER`.
Coordinates and directions share the same row-major layout as
:c:func:`potlibc_evaluate`. The jet buffer has shape
``[order][count][3 * atom_count]``: slice ``jet + (k - 1) * count * 3 *
atom_count`` holds order ``k``. Pass a ``NULL`` :c:struct:`potlibc_jet_report`
when the optional operation count and scale are not needed.

The reported :c:member:`potlibc_jet_report.operations` and
:c:member:`potlibc_jet_report.scale` help estimate roundoff as
``operations * scale * DBL_EPSILON``; they are indicative, not certified bounds.

Always inspect :c:member:`potlibc_metadata.capabilities` before calling an
optional operation. All 387 generated PotLib 2023 adapters in this release
advertise energy, force, and HVP support; the CH5 (``ch52008``) surface also
advertises :c:macro:`POTLIBC_CAP_GRADIENT_JET`.

Errors
------

Every fallible call returns :c:enum:`potlibc_status`.
:c:func:`potlibc_status_string` returns a static diagnostic string suitable
for logs. Failed creation leaves the output handle as ``NULL``. Input and
output arrays may not overlap, and coordinates and directions must be finite.

Ownership and concurrency
-------------------------

A handle owns one private plugin image. Adapters other than CBE restore a
writable-data snapshot before each seed. Different handles may be evaluated
concurrently. Do not overlap two
calls on the same handle, and do not destroy a handle while a call is active.
Create one handle per worker when a surface is shared by a thread pool.
