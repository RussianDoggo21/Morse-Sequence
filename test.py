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

st = SimplexTree(MakeFacesVectorized1(100, 100))

# Créer une instance de MorseSequence
ms = MorseSequence(st)

S = sorted(st.simplices(), key=lambda x: (len(x), x))
F = dict()
for s in S:
    cn = ms.find_node(s)
    F[cn] = 0
S = sorted(st.simplices(), key=lambda s: (F[ms.find_node(s)], len(s)))

print(f" type de S : {type(S)}, type des éléments de S {type(S[0])}")
for s in S : 
    print(f"type de F : {type(F)}, type des clés de F : {type(s)}, type des valeurs de F : {type(F[s])}")
    break
"""
max, n_crit = ms.Max(S, F)

S = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
F = dict()
for s in S:
    F[s] = 0
S = sorted(st.simplices(), key=lambda s: (-F[s], -len(s)))
min, n_crit2 = ms.Min(S, F)


print(f"ms_dec :\n {max},\n n_crit = {n_crit}\n\n")
print(f"ms_crois :\n {min},\n n_crit = {n_crit2}\n\n")
"""
"""
Pour lancer test.py : 
1) Entrer dans l'environnement virtuel : source venv/bin/activate
2) Lancer test.py : python3 test.py
3) Quitter l'environnement virtuel : deactivate
"""