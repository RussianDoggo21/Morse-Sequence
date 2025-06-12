project('morse_sequence', 'cpp',
  version: '0.2',
  default_options: ['cpp_std=c++17']
)

# Find Python installation
py = import('python').find_installation()

# Dependency for pybind11
pybind11_dep = dependency('pybind11', required: true)

# Include directories for simplextree
simplextree_inc = include_directories('simplextree-py/include')

# Source files for simplextree
simplextree_sources = files(
  'simplextree-py/simplextree/_simplextree.cpp'
)

# Define the Python extension module
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
  install: true,
)
