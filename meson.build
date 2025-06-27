project('morse_sequence', 'cpp',
  version : '0.3',
  default_options : ['cpp_std=c++20', 'optimization=3', 'b_lto=true']
)

add_global_arguments(['-march=native', '-flto', '-funroll-loops'], language: 'cpp')

py = import('python').find_installation()          # Python du venv
pybind11_dep = dependency('pybind11', required : true)

# --------------------------------------------------------------------
# chemins utiles
py_site      = run_command(py, '-c', 'import site, sys; print(site.getsitepackages()[0])',
                           check : true).stdout().strip()
pkg_dir      = join_paths(py_site, 'morse_sequence')   # => …/site-packages/morse_sequence
simplextree_inc = include_directories('simplextree-py/include')
tsl_inc = include_directories('external')

# --------------------------------------------------------------------
# sources C++
core_sources = files(
  'src/morse_sequence/_core/bindings.cpp',
  'src/morse_sequence/_core/morse_sequence.cpp',
  'simplextree-py/simplextree/_simplextree.cpp'
)

# --------------------------------------------------------------------
# module natif pybind11 : morse_sequence._core
py.extension_module(
  '_core',        # nom logique du module
  sources            : core_sources,
  include_directories: [
      include_directories('src/morse_sequence/_core'),
      simplextree_inc, 
      tsl_inc
  ],
  dependencies : [pybind11_dep],
  install      : true,
  install_dir  : pkg_dir          # <- on place le .so dans le package
)

# --------------------------------------------------------------------
# installe également __init__.py pour constituer un vrai package
install_data(
  'src/morse_sequence/__init__.py',
  install_dir : pkg_dir
)
