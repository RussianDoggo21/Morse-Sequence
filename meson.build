project('morse_sequence', 'cpp', version : '0.1.0', default_options : ['cpp_std=c++17'])

pybind11_dep = dependency('pybind11')
python = import('python').find_installation(pure: false)

sources = ['morse_sequence.cpp', 'bindings.cpp']

py_ext = python.extension_module('morse_sequence', sources, dependencies : [pybind11_dep], install : true)  # Installe automatiquement le module Python
