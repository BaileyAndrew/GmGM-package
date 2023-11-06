# meson setup build
# meson compile -C build
# WIPE: meson build --wipe

project('GmGM_BaileyA', 'c',
  version : '0.0.1',
  license: 'MIT',
  meson_version: '>=1.1.0',
  default_options : [
    'warning_level=0',
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types -Wno-int-conversion',
    'fortran_args=-Wno-conversion',
  ],
)

add_languages('fortran')

# This shuts up a common irrelevant warning
c_flags = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

message(py.full_path())

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Taken from https://github.com/aburrell/apexpy/blob/main/meson.build
# Unlike distutils, meson doesn't yet include some of the f2py stuff
fortranobject_c = incdir_f2py / 'fortranobject.c'

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: c_flags,
  dependencies: py_dep,
  include_directories: inc_np)

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: inc_np)

# Need to manually run f2py in the fortran directory:
# python -m numpy.f2py fortran_core.f90 -m fortran_core --quiet
py.extension_module('fortran_core',
  [
    'src/GmGM_BaileyA/fortran/fortran_core.f90',
    'src/GmGM_BaileyA/fortran/fortran_core-f2pywrappers2.f90',
    'src/GmGM_BaileyA/fortran/fortran_coremodule.c',
  ],
  c_args : c_flags,
  include_directories: inc_np,
  link_with : fortranobject_lib,
  dependencies : [py_dep, fortranobject_dep],
  install : true,
  subdir : 'GmGM'
)

py.install_sources([
  'src/GmGM_BaileyA/__init__.py',
  'src/GmGM_BaileyA/GmGM.py',
  'src/GmGM_BaileyA/nonparanormal_skeptic.py',
  'src/GmGM_BaileyA/generate_data.py',
  'src/GmGM_BaileyA/regularizers.py'],
  subdir : 'GmGM'
)

message(py.get_install_dir())

# Get build dir
pypathdir = meson.current_build_dir()

message(meson.current_build_dir())

# Check with `meson test -v`
test('extmod',
    py,
    args : files('tests/fortran_imports.py'),
    env : ['PYTHONPATH=' + py.get_install_dir()])