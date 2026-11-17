.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

=======
potlibc
=======

``potlibc`` is a stable C ABI over all 387 potential-energy surface versions
in the PotLib 2023 deposit. Every surface provides energy, analytic force, and
a mixed-forward-dual Hessian-vector product (HVP) without materializing the
Hessian.

Each handle receives a private image of the generated surface plugin. This
isolates legacy Fortran ``COMMON`` and ``SAVE`` state while allowing unrelated
handles to run concurrently without a process-wide evaluation lock.

Coordinates are Bohr, energies Hartree, forces Hartree/Bohr, and HVPs
Hartree/Bohr\ :sup:`2`.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   install
   usage
   rgpot
   accuracy
   provenance

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   surfaces
   design
