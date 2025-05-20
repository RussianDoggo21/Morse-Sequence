from simplextree import SimplexTree
from ..decreasing import morse_seq_decreasing
import numpy as np

"""
st = SimplexTree([[1, 5, 7], [1, 2, 7],    # Top left
                 [2, 7, 9], [2, 3, 9],    # Top middle
                 [3, 5, 9], [1, 3, 5],    # Top right
                 [5, 4, 6], [5, 6, 7],    # Middle left
                 [7, 6, 8], [7, 8, 9],    # Middle center
                 [9, 8, 4], [9, 4, 5],    # Middle right
                 [1, 2, 4], [2, 4, 6],    # Bottom left
                 [2, 3, 6], [3, 6, 8],    # Bottom middle
                 [1, 3, 8], [1, 4, 8]])   # Bottom right
"""
st = SimplexTree([[1,2,3]])
seq, n_crit = morse_seq_decreasing(st)

print(f"n_crit = {n_crit}")
print(f"seq = {seq}") 

# To run the file from the root : python3 -m morse_sequence.tests.test_decreasing