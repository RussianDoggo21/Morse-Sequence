project('morse_sequence', 'cpp', version : '0.1.0', default_options : ['cpp_std=c++20'])

pybind11_dep = dependency('pybind11', required: true)
python = import('python').find_installation(pure: false)

# Ajout explicite du chemin vers les headers Python (nécessaire dans ton cas)
python_inc = include_directories('/usr/include/python3.13')

sources = ['f_sequence.cpp', 'bindings.cpp']

py_ext = python.extension_module(
  'morse_sequence',
  sources,
  dependencies: [pybind11_dep],
  include_directories: python_inc,
  install: true
)
