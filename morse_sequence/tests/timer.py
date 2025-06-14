# File to measure the time complexity of each of the 4 algorithms
# In both pure python and hybrid code

from simplextree import SimplexTree
from morse_sequence import MorseSequence
from ..max import Max
from ..min import Min
from ..increasing import morse_seq_increasing
from ..decreasing import morse_seq_decreasing
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
    list_faces = [10, 20, 50, 60, 75]
    for k in list_faces:
        print(f"\n======================= Cas grille {k} x {k} =======================\n")
        st = SimplexTree(MakeFacesVectorized1(k, k))
        ms = MorseSequence(st)    

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #

        start_dec2 = time.time()
        ms_dec2, n_crit_dec2 = morse_seq_decreasing(st)
        end_dec2 = time.time()
        time_dec2 = end_dec2 - start_dec2
        print(f"décroissante python pur : {time_dec2:.6f} secondes")

        
        start_dec = time.time()
        ms_dec, n_crit_dec = ms.ms_decreasing(st)
        end_dec = time.time()
        time_dec = end_dec - start_dec
        print(f"décroissante C++/python : {time_dec:.6f} secondes")

        print(f"Différence : {time_dec2 - time_dec:.6f} secondes | Ratio (python pur / C++): {time_dec2 / time_dec:.2f}\n")
        
# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #

        start_crois2 = time.time()
        ms_crois2, n_crit_crois2 = morse_seq_increasing(st)
        end_crois2 = time.time()
        time_crois2 = end_crois2 - start_crois2
        print(f"croissante python pur : {time_crois2:.6f} secondes")

        
        start_crois = time.time()
        ms_crois, n_crit_crois = ms.ms_increasing(st)
        end_crois = time.time()
        time_crois = end_crois - start_crois
        print(f"croissante C++/python : {time_crois:.6f} secondes")

        print(f"Différence : {time_crois2 - time_crois:.6f} secondes | Ratio (python pur / C++): {time_crois2 / time_crois:.2f}\n")
        
# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #

        S_max = sorted(st.simplices(), key=lambda x: (len(x), x))
        F_max = {s: 0 for s in S_max}
        S_max = sorted(st.simplices(), key=lambda s: (F_max[s], len(s)))

        start_max2 = time.time()
        max2 = Max(S_max, st, F_max)
        end_max2 = time.time()
        time_max2 = end_max2 - start_max2
        print(f"max python pur : {time_max2:.6f} secondes")

        
        start_max = time.time()
        max, n_crit_max = ms.Max(S_max, F_max)
        end_max = time.time()
        time_max = end_max - start_max
        print(f"max C++/python : {time_max:.6f} secondes")

        print(f"Différence : {time_max2 - time_max:.6f} secondes | Ratio (python pur / C++): {time_max2 / time_max:.2f}\n")
        

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #

        S_min = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
        F_min = {s: 0 for s in S_min}
        S = sorted(st.simplices(), key=lambda s: (-F_min[s], -len(s)))

        start_min2 = time.time()
        min2 = Min(S_min, st, F_min)
        end_min2 = time.time()
        time_min2 = end_min2 - start_min2
        print(f"Min python pur : {time_min2:.6f} secondes")

        
        start_min = time.time()
        min, n_crit_min = ms.Min(S_min, F_min)
        end_min = time.time()
        time_min = end_min - start_min
        print(f"Min C++/python : {time_min:.6f} secondes")

        print(f"Différence : {time_min2 - time_min:.6f} secondes | Ratio (python pur / C++): {time_min2 / time_min:.2f}\n")
        

timer_comparison()

# To run the file from the root : python3 -m morse_sequence.tests.timer

"""
Pour lancer test.py : 
1) Entrer dans l'environnement virtuel : source venv/bin/activate
2) Lancer test.py : python3 main.py
3) Quitter l'environnement virtuel : deactivate
"""