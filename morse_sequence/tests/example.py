from morse_sequence import MorseSequence
from simplextree import SimplexTree

st = SimplexTree([[1,2,3]]) # Creation of simplicial complex via the library SimplexTree
ms = MorseSequence(st) # MorseSequence created

morse_seq_dec, n_crit_dec = ms.ms_decreasing(st) # Computation of a decreasing Morse Sequence on st and its critical simplices
print(f"Critical simplices = {n_crit_dec}") # Critical simplices = 1
print(f"Decreasing Morse Sequence = {morse_seq_dec}") # Decreasing Morse Sequence = [((2, 3), (1, 2, 3)), ((3,), (1, 3)), ((1,), (1, 2)), [(2,)]]

morse_seq_inc, n_crit_inc = ms.ms_increasing(st) # Computation of an increasing Morse Sequence on st and its critical simplices
print(f"Critical simplices = {n_crit_inc}") # Critical simplices = 1
print(f"Increasing Morse Sequence = {morse_seq_inc}") # Increasing Morse Sequence = [[1], ([3], [1, 3]), ([2], [2, 3]), ([1, 2], [1, 2, 3])]