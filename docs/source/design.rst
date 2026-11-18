.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

Generated adapter design
========================

Inventory and generation
------------------------

``tools/surface_inventory.py`` classifies all 403 deposited Fortran files and
selects 387 surface roots. At Meson configure time, that inventory determines
one custom generation target and one shared plugin per surface.

The generator chooses an adapter from the declared PotLib interface family or
the parsed callable signature. Parameter files are embedded so plugin
instances do not depend on a mutable working directory. The generated plugin
exports a small internal ABI consumed only by ``libpotlibc``; applications use
the stable public ABI in ``potlibc.h``.

Dual propagation
----------------

Each transformed scalar carries primal, coordinate, direction, and mixed
lanes. The coordinate lane supplies analytic gradients, while the mixed lane
supplies the directional derivative of that gradient. This computes
:math:`H(x)v` with storage linear in the coordinate count.

State isolation
---------------

On handle creation, ``libpotlibc`` copies the selected plugin into an
anonymous file and loads that unique image. Initialization runs once.
Generated adapters other than CBE (``ch52008``) then snapshot writable ELF
segments, and each coordinate or direction seed restores that snapshot
before invoking transformed surface code. CBE keeps the deposited analytic
``pot`` path and does not restore a snapshot. Its isolation is still the
private plugin image.

Symbolic references are bound within each plugin. Distinct handles therefore
do not share the surface's Fortran ``COMMON`` blocks or ``SAVE`` variables,
while standard runtime libraries remain shared in the base namespace.

Platform boundary
-----------------

The current loader depends on Linux ``memfd_create``, ``/proc/self/fd``, and
the thread-safe glibc ``dlopen``/``dlclose`` implementation. Other platforms
are not supported by this release, even though the public header contains
portable symbol-visibility declarations.
