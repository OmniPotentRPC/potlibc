.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

C API reference
===============

The public ABI is declared entirely in ``include/potlibc.h``. This reference
is generated from that header by Doxygen and Breathe.

Types
-----

.. doxygentypedef:: potlibc_handle
   :project: potlibc

.. doxygenenum:: potlibc_status
   :project: potlibc

.. doxygenenum:: potlibc_capability
   :project: potlibc

.. doxygenstruct:: potlibc_metadata
   :project: potlibc
   :members:

Functions
---------

.. doxygenfunction:: potlibc_create
   :project: potlibc

.. doxygenfunction:: potlibc_create_sized
   :project: potlibc

.. doxygenfunction:: potlibc_destroy
   :project: potlibc

.. doxygenfunction:: potlibc_get_metadata
   :project: potlibc

.. doxygenfunction:: potlibc_evaluate
   :project: potlibc

.. doxygenfunction:: potlibc_energy_directional_derivative
   :project: potlibc

.. doxygenfunction:: potlibc_hessian_vector_product
   :project: potlibc

.. doxygenfunction:: potlibc_status_string
   :project: potlibc
