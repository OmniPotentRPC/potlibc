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

Corrected gradients
-------------------

Gradients are consistent with the returned energy by default. The deposited
CH5 2008 CBE source uses a rounded gradient conversion, so potlibc applies its
cataloged correction to both forces and Hessian-vector products. Call
:c:func:`potlibc_set_gradient_mode` with
:c:enumerator:`POTLIBC_GRADIENT_RAW_UPSTREAM` only when reproducing legacy
upstream gradient bits is more important than derivative consistency.

Hessian-vector products
-----------------------

:c:func:`potlibc_hessian_vector_product` computes
:math:`H(x)v = \nabla^2 E(x)v` directly with mixed forward dual numbers. It
does not allocate or expose a dense Hessian.

.. literalinclude:: ../../examples/hvp.c
   :language: c
   :linenos:

Always inspect :c:member:`potlibc_metadata.capabilities` before calling an
optional operation. All 387 generated PotLib 2023 adapters in this release
advertise energy, force, and HVP support.

Errors
------

Every fallible call returns :c:enum:`potlibc_status`.
:c:func:`potlibc_status_string` returns a static diagnostic string suitable
for logs. Failed creation leaves the output handle as ``NULL``. Input and
output arrays may not overlap, and coordinates and directions must be finite.

Ownership and concurrency
-------------------------

A handle owns one private plugin image and its restored writable-state
snapshot. Different handles may be evaluated concurrently. Do not overlap two
calls on the same handle, and do not destroy a handle while a call is active.
Create one handle per worker when a surface is shared by a thread pool.
