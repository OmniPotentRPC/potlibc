# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause
project(
  'potlibc',
  ['c', 'fortran'],
  version: '0.2.0',
  meson_version: '>=1.1.0',
  default_options: ['c_std=c11', 'warning_level=3'],
)

python = import('python').find_installation('python3')
pkg = import('pkgconfig')
cc = meson.get_compiler('c')
lapack = dependency('lapack')
threads = dependency('threads')
dl = cc.find_library('dl')
catalog = files('tools/catalog.py')
inventory = files('tools/surface_inventory.py')
generator = files('tools/generate_surface.py')
surface_metadata = files('data/surfaces.json')
vendor = meson.current_source_dir() / 'upstream/potlib2023'
plugin_install_dir = get_option('libdir') / 'potlibc'

inventory_result = run_command(
  python,
  inventory,
  vendor,
  meson.current_build_dir() / 'surface-inventory.json',
  '--meson',
  check: true,
)
surface_rows = inventory_result.stdout().strip().split('\n')
surface_modules = []
surface_ids = []
surface_test_cases = []
foreach row : surface_rows
  columns = row.split('|')
  surface_id = columns[0]
  surface_ids += [surface_id]
  surface_test_cases += [surface_id + '=' + columns[2]]
  vendor_source = files('upstream/potlib2023' / columns[1])
  generated_outputs = [
    surface_id + '-surface.f90',
    surface_id + '-utility.f90',
    surface_id + '-adapter.f90',
    surface_id + '-data.c',
    surface_id + '-data.bin',
  ]
  generated_command = [
    python, generator, vendor, surface_id,
    '@OUTPUT0@', '@OUTPUT1@', '@OUTPUT2@', '@OUTPUT3@', '@OUTPUT4@',
      '--metadata-catalog', surface_metadata,
  ]
  generator_inputs = [
    generator,
    inventory,
    catalog,
    surface_metadata,
    vendor_source,
    files('upstream/potlib2023/utilities/utility.f'),
  ]
  if surface_id == 'n4pes-l-imls-g2-v1'
    generator_inputs += files(
      'upstream/potlib2023/src/N4/n4pes_L-IMLS-G2_v1.tar',
    )
  endif
  generated = custom_target(
    'generate-' + surface_id,
    input: generator_inputs,
    output: generated_outputs,
    command: generated_command,
  )
  module = shared_module(
    'potlibc-' + surface_id,
    files('src/dual2_mod.f90', 'src/plugin_state.c') + [
      generated[0],
      generated[1],
      generated[2],
      generated[3],
    ],
    fortran_args: (
      surface_id == 'ch52008' ? ['-O3'] : []
    ) + (
      surface_id in [
        'phsch3-ia', 'phsch3-s-1-8', 'phsch3-s-2-0', 'phsch3-ii',
        'phsch3-s-2-5', 'phsch3-iii', 'phsch3-s-2-2',
        'phsch3-aprp', 'ch4t8',
      ] ? ['-O0'] : []
    ) + (
      surface_id in ['h3dmbe', 'sioh', 'oh2m3a2001'] ?
        ['-finit-derived', '-finit-real=zero'] : []
    ) + (
      surface_id in [
        'phsch3-ia', 'phsch3-s-1-8', 'phsch3-s-2-0', 'phsch3-ii',
        'phsch3-s-2-5', 'phsch3-iii', 'ch3nh2', 'phsch3-s-2-2',
        'phsch3-aprp', 'phsh-aprp', 'ch4t8',
      ] ?
        ['-finit-derived', '-finit-real=zero'] : []
    ) + (
      surface_id.startswith('c2h6') ?
        ['-finit-derived', '-finit-real=zero'] : []
    ) + [
      '-std=legacy',
      '-ffree-line-length-none',
      '-fno-implicit-none',
      '-fallow-argument-mismatch',
      '-fno-align-commons',
      '-fPIC',
    ],
    link_args: ['-Wl,-z,defs', '-Wl,-Bsymbolic'],
    dependencies: [lapack, dl],
    override_options: ['warning_level=0'],
    install: true,
    install_dir: plugin_install_dir,
    install_rpath: get_option('libdir'),
  )
  surface_modules += module
endforeach

potlibc_inc = include_directories('include')
potlibc = library(
  'potlibc',
  'src/potlibc.c',
  include_directories: potlibc_inc,
  c_args: [
    '-DPOTLIBC_BUILDING_LIBRARY',
    '-DPOTLIBC_PLUGIN_DIR="@0@"'.format(
      get_option('prefix') / plugin_install_dir
    ),
  ],
  dependencies: [dl, threads],
  gnu_symbol_visibility: 'hidden',
  install: true,
  # Compatible plugin and C ABI addition. SONAME stays libpotlibc.so.0.
  version: meson.project_version(),
  soversion: '0',
)

install_headers('include/potlibc.h')
pkg.generate(
  potlibc,
  description: 'Isolated C API for POTLIB potential-energy surfaces',
  filebase: 'potlibc',
  name: 'potlibc',
  subdirs: '.',
)

if get_option('examples')
  subdir('examples')
endif
subdir('tests')
