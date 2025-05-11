from simplextree import SimplexTree
from morse_sequence import MorseSequence
import numpy as np
import time


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

def timer_comparison():
    list_faces = [10, 20, 50, 75, 100]
    for k in list_faces:
        st = SimplexTree(MakeFacesVectorized1(k, k))
        ms = MorseSequence(st)    

        start_dec = time.time()
        ms_dec, n_crit_dec = ms.morse_seq_decrois(st)
        end_dec = time.time()
        print(f"Temps d'exécution décroissante, st = MakeFacesVectorized1({k}, {k}): {end_dec - start_dec:.6f} secondes\n\n")

        start_crois = time.time()
        ms_crois, n_crit_crois = ms.morse_seq_crois(st)
        end_crois = time.time()
        print(f"Temps d'exécution croissante, st = MakeFacesVectorized1({k}, {k}): {end_crois - start_crois:.6f} secondes\n\n")

        S_max = sorted(st.simplices(), key=lambda x: (len(x), x))
        F_max = {s: 0 for s in S_max}
        S_max = sorted(st.simplices(), key=lambda s: (F_max[s], len(s)))
        start_max = time.time()
        max, n_crit_max = ms.Max(S_max, F_max)
        end_max = time.time()
        print(f"Temps d'exécution max, st = MakeFacesVectorized1({k}, {k}): {end_max - start_max:.6f} secondes\n\n")

        S_min = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
        F_min = {s: 0 for s in S_min}
        S = sorted(st.simplices(), key=lambda s: (-F_min[s], -len(s)))
        start_min = time.time()
        min, n_crit_min = ms.Min(S_min, F_min)
        end_min = time.time()
        print(f"Temps d'exécution min, st = MakeFacesVectorized1({k}, {k}): {end_min - start_min:.6f} secondes\n\n")    

timer_comparison()

""""
st = SimplexTree(MakeFacesVectorized1(10, 10))
#st = SimplexTree([[1,2,3]])
#st = SimplexTree([[1, 5, 7], [1, 2, 7], [2, 7, 9], [2, 3, 9], [3, 5, 9], [1, 3, 5], [5, 4, 6], [5, 6, 7], [7, 6, 8],
                   #[7, 8, 9], [9, 8, 4], [9, 4, 5], [1, 2, 4], [2, 4, 6], [2, 3, 6], [3, 6, 8], [1, 3, 8], [1, 4, 8]])

# Créer une instance de MorseSequence
ms = MorseSequence(st)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
start_dec = time.time()
ms_dec, n_crit_dec = ms.morse_seq_decrois(st)
end_dec = time.time()
#print(f"decroissante :\n{ms_dec},\n n_crit = {n_crit_dec}\n")
print(f"Temps d'exécution décroissante : {end_dec - start_dec:.6f} secondes\n\n")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
start_crois = time.time()
ms_crois, n_crit_crois = ms.morse_seq_crois(st)
end_crois = time.time()
#print(f"croissante :\n{ms_crois},\n n_crit = {n_crit_crois}\n")
print(f"Temps d'exécution croissante : {end_crois - start_crois:.6f} secondes\n\n")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
S_max = sorted(st.simplices(), key=lambda x: (len(x), x))
# Dictionnaire F avec les pointeurs comme clés
F_max = {s: 0 for s in S_max}
# Tri selon F et la dimension du simplexe (exposée en C++)
S_max = sorted(st.simplices(), key=lambda s: (F_max[s], len(s)))

start_max = time.time()
max, n_crit_max = ms.Max(S_max, F_max)
end_max = time.time()
#print(f"max :\n{max},\n n_crit = {n_crit_max}\n")
print(f"Temps d'exécution max: {end_max - start_max:.6f} secondes\n\n")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #

S_min = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
F_min = {s: 0 for s in S_min}
S = sorted(st.simplices(), key=lambda s: (-F_min[s], -len(s)))

start_min = time.time()
min, n_crit_min = ms.Min(S_min, F_min)
end_min = time.time()
#print(f"min :\n{min},\n n_crit = {n_crit_min}\n")
print(f"Temps d'exécution min: {end_min - start_min:.6f} secondes\n\n")
"""
"""
Pour lancer test.py : 
1) Entrer dans l'environnement virtuel : source venv/bin/activate
2) Lancer test.py : python3 main.py
3) Quitter l'environnement virtuel : deactivate
"""