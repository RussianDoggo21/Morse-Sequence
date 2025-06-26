from simplextree import SimplexTree
from morse_sequence import MorseSequence
from morse_sequence._core import SimplexBatch

from ..increasing import morse_seq_increasing
import numpy as np


def MakeFacesVectorized1(Nr,Nc):

    out = np.empty((Nr-1,Nc-1,2,3),dtype=int)

    r = np.arange(Nr*Nc).reshape(Nr,Nc)

    out[:,:, 0,0] = r[:-1,:-1]
    out[:,:, 1,0] = r[:-1,1:]
    out[:,:, 0,1] = r[:-1,1:]

    out[:,:, 1,1] = r[1:,1:]
    out[:,:, :,2] = r[1:,:-1,None]

    out.shape =(-1,3)
    return out




st = SimplexTree([[1, 5, 7], [1, 2, 7],    # Top left
                 [2, 7, 9], [2, 3, 9],    # Top middle
                 [3, 5, 9], [1, 3, 5],    # Top right
                 [5, 4, 6], [5, 6, 7],    # Middle left
                 [7, 6, 8], [7, 8, 9],    # Middle center
                 [9, 8, 4], [9, 4, 5],    # Middle right
                 [1, 2, 4], [2, 4, 6],    # Bottom left
                 [2, 3, 6], [3, 6, 8],    # Bottom middle
                 [1, 3, 8], [1, 4, 8]])   # Bottom right




seq, n_crit = morse_seq_increasing(st)

print(f"n_crit = {n_crit}")
print(f"seq : {seq}\n")
"""
st = SimplexTree([[1,2,3]])
ms = MorseSequence(st)    

simplices = st.simplices()

# batch S (cosimplicial complex)
batch_S = SimplexBatch.from_python(simplices, st) # no weights

# batch F (dictionnary simplex : weight)
simp_with_weights = [(list(s), 0) for s in simplices]     # [(simplexe, weight), ...]
batch_F = SimplexBatch.from_python(simp_with_weights, st) # weights.size() == nodes.size()

max, n_crit_max = ms.Max(batch_S, batch_F)



"""
# To run the file from the root : python3 -m src.python.python_tests.test_increasing