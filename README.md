# Morse-Sequence

Repository created for the report of an internship as an assistant of research at ESIEE PARIS under Mr. NAJMAN L.
Based on the library SimplexTree to compute Morse Sequences and F-Sequences from a simplicial complex.

## Prerequisites

You will first need to install the library SimplexTree either via pip or locally.
Refer to the following for more information: [simplextree-py](https://github.com/peekxc/simplextree-py)

C++ and Python compilers are also necessary.

## How to install morse_sequence

### Clone with submodules

```bash
git clone --recurse-submodules https://github.com/ton-compte/morse_frame.git
cd morse_frame
```

### Installation (standard)

To build and install the package, run:

```bash
pip install .
```

This will compile the C++ code via Meson and install the package.

### Installation in development mode (editable)

If you want to install the package in editable mode to develop locally:

```bash
pip install --no-build-isolation --editable .
```

This allows you to modify the code without reinstalling after every change.

## Usage example

```python
from morse_sequence import MorseSequence

# Your code here
```
