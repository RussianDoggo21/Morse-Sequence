from simplextree import SimplexTree
from morse_sequence import MorseSequence
import numpy as np

# Créer une instance de SimplexTree
# tree = SimplexTree([[1, 5, 7], [1, 2, 7], [2, 7, 9], [2, 3, 9], [3, 5, 9], [1, 3, 5], [5, 4, 6], [5, 6, 7], [7, 6, 8],
    #               [7, 8, 9], [9, 8, 4], [9, 4, 5], [1, 2, 4], [2, 4, 6], [2, 3, 6], [3, 6, 8], [1, 3, 8], [1, 4, 8]])


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

tree = SimplexTree(MakeFacesVectorized1(100, 100))

# Créer une instance de MorseSequence
ms = MorseSequence(tree)

# Appeler les méthodes
ms_dec, n_crit = ms.morse_seq_dec(tree)
ms_crois, n_crit2 = ms.morse_seq_crois(tree)


print(f"ms_dec :\n {ms_dec},\n n_crit = {n_crit}\n\n")
print(f"ms_crois :\n {ms_crois},\n n_crit = {n_crit2}\n\n")

"""
Pour lancer test.py : 
1) Entrer dans l'environnement virtuel : source venv/bin/activate
2) Lancer test.py : python3 test.py
3) Quitter l'environnement virtuel : deactivate
"""