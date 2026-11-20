.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

Provenance and licensing
========================

Vendor source
-------------

``upstream/potlib2023`` retains the PotLib 2023 deposit:

  Y. Shu, Z. Varga, A. Jasper, J. Espinosa-García, J. C. Corchado, and
  D. G. Truhlar, "PotLib 2023: New version of a potential energy surface
  library for chemical systems," Mendeley Data, version 1 (2023),
  `doi:10.17632/chxxhnnt4p.1 <https://doi.org/10.17632/chxxhnnt4p.1>`_.

The repository documents local numerical corrections to vendored files in
``PROVENANCE.md``, which records a SHA-256 digest for every vendored file.
``data/surfaces.json`` maps all 403 deposited Fortran files to 387 callable
surface roots and 16 owned support files.

Build-time transforms
---------------------

The build reads the vendor tree and writes derived files only into the Meson
build directory. For each surface,
``tools/generate_surface.py`` emits:

* a normalized generated Fortran surface and utility unit;
* a family-specific C adapter;
* embedded parameter data; and
* a private surface plugin.

The transform normalizes legacy fixed-form declarations, applies documented
modern-compiler repairs, replaces supported arithmetic with mixed forward
duals, and embeds packaged parameter records. A small, enumerated set of
legacy implementations receives compiler initialization or optimization
flags in ``meson.build``. These operations do not rewrite the deposited
files.

Licences
--------

The complete ``upstream/potlib2023`` vendor tree is distributed verbatim
under Apache-2.0. The exact text is in ``LICENSES/Apache-2.0.txt`` and the
file-level declaration is in ``.reuse/dep5``.

New ``potlibc`` source code, generated catalog metadata, examples, and
documentation are BSD-3-Clause; see ``LICENSES/BSD-3-Clause.txt``. Generated
build products are derived from both the Apache-2.0 vendor inputs and the
BSD-3-Clause transform and adapter code. Downstream distributors must retain
the applicable notices and licence texts.
