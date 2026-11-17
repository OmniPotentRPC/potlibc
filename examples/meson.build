# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause

example_env = environment()
example_env.set('POTLIBC_PLUGIN_DIR', meson.current_build_dir() / '..')

foreach example_name : ['energy', 'hvp']
  example = executable(
    'potlibc-' + example_name,
    example_name + '.c',
    include_directories: potlibc_inc,
    link_with: potlibc,
    install: false,
  )
  test('example ' + example_name, example, env: example_env)
endforeach
