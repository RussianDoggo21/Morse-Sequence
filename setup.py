from setuptools import setup, find_packages

setup(
    name='morse_sequence',
    version='0.2',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=['pybind11'],
)
