from simplextree import SimplexTree
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
seq, n_crit = morse_seq_increasing(st)

print(n_crit)
print(seq, "\n")


# To run the file from the root : python3 -m morse_sequence.tests.test_increasing