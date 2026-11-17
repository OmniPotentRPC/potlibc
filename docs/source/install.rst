.. SPDX-FileCopyrightText: 2026 Rohit Goswami
.. SPDX-License-Identifier: BSD-3-Clause

Install and build
=================

Requirements
------------

``potlibc`` currently targets Linux because per-handle isolation uses
``memfd_create`` and ``/proc/self/fd``. A build needs:

* Meson 1.1 or newer and Ninja;
* a C11 compiler and a GNU-compatible Fortran compiler;
* LAPACK, POSIX threads, and ``libdl`` development files; and
* Python 3 for the catalog and build-time source transforms.

On Debian or Ubuntu, the corresponding packages are:

.. code-block:: console

   $ sudo apt-get install gcc gfortran liblapack-dev meson ninja-build python3

Build and test
--------------

.. code-block:: console

   $ meson setup build --buildtype=release
   $ meson compile -C build
   $ meson test -C build --print-errorlogs

The initial compile generates and builds 387 isolated surface plugins. Use
``meson setup --wipe build`` after changing compilers or LAPACK providers.
ThreadSanitizer is exposed as a separate target when the compiler supports it:

.. code-block:: console

   $ meson compile -C build test-tsan

Install
-------

.. code-block:: console

   $ meson setup build --prefix="$HOME/.local" --buildtype=release
   $ meson compile -C build
   $ meson install -C build

This installs ``libpotlibc``, ``potlibc.h``, the ``potlibc.pc`` pkg-config
file, and one plugin per catalog surface under the configured library
directory. Consumers can use:

.. code-block:: console

   $ cc example.c $(pkg-config --cflags --libs potlibc)

At runtime, an installed library uses its configured plugin directory.
``POTLIBC_PLUGIN_DIR`` may override that location for an uninstalled build or
packaging layout.

Build the documentation
-----------------------

Install Doxygen plus the Python packages in ``docs/requirements.txt``, then:

.. code-block:: console

   $ python3 -m pip install -r docs/requirements.txt
   $ python3 docs/build.py

The HTML output is written to ``docs/_build/html``. Doxygen and Sphinx warnings
are errors.
