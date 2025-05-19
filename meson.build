project('morse_sequence', 'cpp',
  version: '0.1',
  default_options: ['cpp_std=c++17']
)

py = import('python').find_installation()

pybind11_dep = dependency('pybind11', required: true)

# Include simplextree headers
simplextree_inc = include_directories('simplextree-py/include')

# Sources de simplextree
simplextree_sources = files(
  'simplextree-py/simplextree/_simplextree.cpp',
  # Ajoute d'autres fichiers .cpp si nécessaire
)

# Extension Python
py.extension_module('morse_sequence',
  sources: [
    'src/morse_sequence/bindings.cpp',
    'src/morse_sequence/morse_sequence.cpp',
    simplextree_sources
  ],
  include_directories: [
    include_directories('src/morse_sequence'),
    simplextree_inc
  ],
  dependencies: [pybind11_dep],
  install: true
)