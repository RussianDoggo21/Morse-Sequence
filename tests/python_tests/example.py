import morse_sequence._core
print(morse_sequence._core.__file__)
print(dir(morse_sequence._core))
print(dir(morse_sequence._core.MorseSequence))

"""
#from morse_sequence import MorseSequence
import morse_sequence._core as _core
from simplextree import SimplexTree

st = SimplexTree([[1,2,3]]) 
ms = _core.MorseSequence(st) # MorseSequence created

morse_seq_dec, n_crit_dec = ms.decreasing() # Computation of a decreasing Morse Sequence on st and its critical simplices
print(f"Critical simplices = {n_crit_dec}") # Critical simplices = 1
print(f"Decreasing Morse Sequence = {morse_seq_dec}") # Decreasing Morse Sequence = [((2, 3), (1, 2, 3)), ((3,), (1, 3)), ((1,), (1, 2)), [(2,)]]

morse_seq_inc, n_crit_inc = ms.increasing() # Computation of an increasing Morse Sequence on st and its critical simplices
print(f"Critical simplices = {n_crit_inc}") # Critical simplices = 1
print(f"Increasing Morse Sequence = {morse_seq_inc}") # Increasing Morse Sequence = [[1], ([3], [1, 3]), ([2], [2, 3]), ([1, 2], [1, 2, 3])]

# To run the file from the root : python3 -m tests.python_tests.example
"""