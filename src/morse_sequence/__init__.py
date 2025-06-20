"""
Public API of the 'morse_sequence' package
The native module '_core' implements the 'MorseSequence' class
"""

# Charge tout depuis l'extension C++
from ._core import *          # exporte MorseSequence, reference_map, ...

# Définis __all__ pour rendre explicites les symboles publics (optionnel)
__all__ = ['MorseSequence', 'reference_map', 'coreference_map']
