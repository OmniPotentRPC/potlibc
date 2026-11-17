# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause

threads = dependency('threads')
test_env = environment()
test_env.set('POTLIBC_PLUGIN_DIR', meson.current_build_dir() / '..')
test_env.set(
  'TSAN_OPTIONS',
  'suppressions=' + (meson.current_source_dir() / 'tsan.supp'),
)

dual_edge_test = executable(
  'test-dual2',
  ['test_dual2.f90', '../src/dual2_mod.f90'],
  dependencies: lapack,
)
test('dual intrinsic edge cases', dual_edge_test)

c_api_test = executable(
  'test-c-api',
  'test_c_api.c',
  include_directories: potlibc_inc,
  link_with: potlibc,
  dependencies: cc.find_library('m', required: false),
)
test('C API and CBE exact parity', c_api_test, env: test_env)

fortran_api_test = executable(
  'test-fortran-api',
  'test_fortran_api.f90',
  include_directories: potlibc_inc,
  link_with: potlibc,
)
test('Fortran API and CBE exact parity', fortran_api_test, env: test_env)

thread_test = executable(
  'test-threads',
  'test_threads.c',
  include_directories: potlibc_inc,
  link_with: potlibc,
  dependencies: threads,
)
test('parallel isolated handles', thread_test, env: test_env)

overhead_test = executable(
  'test-overhead',
  'test_overhead.c',
  include_directories: potlibc_inc,
  link_with: potlibc,
  dependencies: cc.find_library('dl'),
)
test('per-evaluation overhead', overhead_test, env: test_env)

consistency_test = executable(
  'test-consistency',
  'test_consistency.c',
  include_directories: potlibc_inc,
  link_with: potlibc,
  dependencies: cc.find_library('m'),
)
foreach surface_case : surface_test_cases
  test(
    'dual consistency ' + surface_case,
    consistency_test,
    args: [surface_case],
    env: test_env,
  )
endforeach

test(
  'generated catalog consistency',
  python,
  args: [
    files('test_catalog.py'),
    inventory,
    vendor,
    files('../data/surfaces.json'),
    meson.current_build_dir() / '..',
  ],
)

if cc.has_argument('-fsanitize=thread') and cc.has_link_argument('-fsanitize=thread')
  tsan_test = executable(
    'test-threads-tsan',
    ['test_threads.c', '../src/potlibc.c'],
    include_directories: potlibc_inc,
    c_args: [
      '-fsanitize=thread',
      '-DPOTLIBC_BUILDING_LIBRARY',
      '-DPOTLIBC_PLUGIN_DIR="@0@"'.format(meson.current_build_dir() / '..'),
    ],
    link_args: ['-fsanitize=thread'],
    dependencies: [threads, cc.find_library('dl')],
  )
  run_target(
    'test-tsan',
    command: tsan_test,
    env: test_env,
  )
endif
