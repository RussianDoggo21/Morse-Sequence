# File to measure the time complexity of each of the 4 algorithms
# In both pure python and hybrid code

from simplextree import SimplexTree
from morse_sequence import MorseSequence
#from morse_sequence._core import SimplexBatch, batch_to_py_simplices

from ..max import Max
from ..min import Min
from ..increasing import morse_seq_increasing
from ..decreasing import morse_seq_decreasing
import numpy as np
import time
from collections import defaultdict

#import cProfile, pstats, io
#from pstats import SortKey

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

""""
def pad_and_pack_S(S_list, dim):

    #Convert list of simplexes (lists of ints) all of length dim+1
    #into a numpy 2D array (shape: (len(S_list), dim+1))
    # S_list is list of tuples/lists each with length = dim+1
    return np.array(S_list, dtype=np.int32)

def pad_and_pack_F(F_list, dim):

    #Convert list of (simplex, value) pairs,
    #where simplex length = dim+1 and value is int,
    #into numpy 2D array with last column = value
    #shape = (len(F_list), dim+2)

    packed = []
    for sigma, val in F_list:
        # sigma is a tuple/list of vertices, val is int
        packed.append(list(sigma) + [val])
    return np.array(packed, dtype=np.int32)
"""

def pad_and_pack_S_general(S_list):
    """
    Pack a list of simplexes (of varying lengths) into a 2D np.array,
    padding with -1 to make rows uniform.
    """
    max_dim = max(len(s) for s in S_list)
    padded = [list(s) + [-1] * (max_dim - len(s)) for s in S_list]
    return np.array(padded, dtype=np.int32)



def pad_and_pack_F_general(F_dict):
    """
    Convert dict of {simplex tuple: int} into a 2D array
    with padded simplexes + value as last column.
    """
    max_dim = max(len(s) for s in F_dict)
    return np.array([
        list(s) + [-1] * (max_dim - len(s)) + [val]
        for s, val in F_dict.items()
    ], dtype=np.int32)


def check_inputs(S, F):
    print("\n[check_inputs]")
    for d in S:
        print(f"  S[{d}] shape: {S[d].shape}, dtype: {S[d].dtype}")
        if len(S[d].shape) != 2:
            raise ValueError(f"S[{d}] doit être 2D")
        if S[d].shape[1] < 1:
            raise ValueError(f"S[{d}] a des lignes trop courtes")

    for d in F:
        print(f"  F[{d}] shape: {F[d].shape}, dtype: {F[d].dtype}")
        if len(F[d].shape) != 2:
            raise ValueError(f"F[{d}] doit être 2D")
        if F[d].shape[1] < 2:
            raise ValueError(f"F[{d}] doit avoir au moins une valeur + un poids (shape[1] >= 2)")
        if d not in S:
            print(f"  ⚠️ Attention : F[{d}] n’a pas de S[{d}] correspondant.")

    if S.keys() != F.keys():
        print("  ⚠️ Clés S ≠ Clés F : cela peut poser problème !")





def timer_comparison():
    list_faces = [10, 20, 50, 60, 75, 100]
    for k in list_faces:
        print(f"\n======================= Cas grille {k} x {k} =======================\n")
        st = SimplexTree(MakeFacesVectorized1(k, k))
        ms = MorseSequence(st)    

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        """
        start_dec2 = time.time()
        ms_dec2, n_crit_dec2 = morse_seq_decreasing(st)
        end_dec2 = time.time()
        time_dec2 = end_dec2 - start_dec2
        print(f"décroissante python pur : {time_dec2:.6f} secondes")

        
        start_dec = time.time()
        ms_dec, n_crit_dec = ms.decreasing(st)
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
        ms_crois, n_crit_crois = ms.increasing(st)
        end_crois = time.time()
        time_crois = end_crois - start_crois
        print(f"croissante C++/python : {time_crois:.6f} secondes")

        print(f"Différence : {time_crois2 - time_crois:.6f} secondes | Ratio (python pur / C++): {time_crois2 / time_crois:.2f}\n")
        """
# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        S_max1 = sorted(st.simplices(), key=lambda x: (len(x), x))
        F_max1 = {s: 0 for s in S_max1}
        S_max1 = sorted(st.simplices(), key=lambda s: (F_max1[s], len(s)))


        start_max2 = time.time()
        max2 = Max(S_max1, st, F_max1)
        end_max2 = time.time()
        time_max2 = end_max2 - start_max2
        print(f"max python pur : {time_max2:.6f} secondes")

        """
        simplices = st.simplices()

        # batch S (cosimplicial complex)
        batch_S = SimplexBatch.from_python(simplices, st) # no weights

        # batch F (dictionnary simplex : weight)
        simp_with_weights = [(list(s), 0) for s in simplices]     # [(simplexe, weight), ...]
        batch_F = SimplexBatch.from_python(simp_with_weights, st) # weights.size() == nodes.size()
        """

        """
        # Convert all simplices to list[int] instead of tuple 
        simplices = [list(s) for s in st.simplices()]
        # Sort in increasing dimension
        S_max2 = sorted(simplices, key=lambda x: (len(x), x))
        # Build F_max as a list of pairs (simplex, weight)
        F_max2 = [[s, 0] for s in S_max2]  # F_max: list[list[int], int]
        # Rebuild a dictionary temporarily for sorting
        F_max_dict2 = {tuple(s): w for s, w in F_max2}
        # Re-sort S_max based on F value and dimension
        S_max2 = sorted(S_max2, key=lambda s: (F_max_dict2[tuple(s)], len(s)))
        """

        """
        # Group S and F by dimension
        S_by_dim_max = defaultdict(list)
        F_by_dim_max = defaultdict(list)

        for sigma in S_max1:
                S_by_dim_max[len(sigma)-1].append(sigma)
        for sigma, val in F_max1.items():
                F_by_dim_max[len(sigma)-1].append((sigma, val))

        # Prepare numpy arrays for each dimension
        S_arrays_max = {}
        F_arrays_max = {}

        for dim, simplices in S_by_dim_max.items():
                S_arrays_max[dim] = pad_and_pack_S(simplices, dim)

        for dim, simplices_vals in F_by_dim_max.items():
                F_arrays_max[dim] = pad_and_pack_F(simplices_vals, dim)

        # Then pass S_arrays_max[dim] and F_arrays_max[dim] to C++ binding functions
        # for each dimension dim.

        for dim in S_arrays_max:
                assert isinstance(S_arrays_max[dim], np.ndarray), f"S_arrays_max[{dim}] is not an ndarray"
                assert isinstance(F_arrays_max[dim], np.ndarray), f"F_arrays_max[{dim}] is not an ndarray"
                assert S_arrays_max[dim].dtype == np.int32, f"S_arrays_max[{dim}] has dtype {S_arrays_max[dim].dtype}"
                assert F_arrays_max[dim].dtype == np.int32, f"F_arrays_max[{dim}] has dtype {F_arrays_max[dim].dtype}"
                assert S_arrays_max[dim].ndim == 2, f"S_arrays_max[{dim}] is not 2D"
                assert F_arrays_max[dim].ndim == 2, f"F_arrays_max[{dim}] is not 2D"
        check_inputs(S_arrays_max, F_arrays_max)
        """

        S_buffer_max = pad_and_pack_S_general(S_max1)
        F_buffer_max = pad_and_pack_F_general(F_max1)

        start_max = time.time()
        max, n_crit_max = ms.Max(S_buffer_max, F_buffer_max)
        end_max = time.time()
        time_max = end_max - start_max
        print(f"max C++/python : {time_max:.6f} secondes")

        print(f"Différence : {time_max2 - time_max:.6f} secondes | Ratio (python pur / C++): {time_max2 / time_max:.2f}\n")
        

# ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        
        S_min1 = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
        F_min1 = {s: 0 for s in S_min1}
        S_min1 = sorted(st.simplices(), key=lambda s: (-F_min1[s], -len(s)))

        start_min2 = time.time()
        min2 = Min(S_min1, st, F_min1)
        end_min2 = time.time()
        time_min2 = end_min2 - start_min2
        print(f"Min python pur : {time_min2:.6f} secondes")

        """
        simplices = st.simplices()                   

        simp_w  = [(list(s), 0) for s in simplices]
        batch_F = SimplexBatch.from_python(simp_w, st)

        batch_S = SimplexBatch.from_python(simplices, st)
        # Sort on the C++ simplices (quick)
        order_S = sorted(range(len(batch_S.simplices)),key=lambda i: (-len(batch_S.simplices[i]), batch_S.simplices[i]))
        # ↪ conversion C++ → Python 
        S_min_py = batch_to_py_simplices(batch_S, st)
        S_min_py_sorted = [S_min_py[i] for i in order_S]   
        batch_S = SimplexBatch.from_python(S_min_py_sorted, st)
        """

        """
        # Convert all simplices to list[int] instead of tuple 
        simplices = [list(s) for s in st.simplices()]
        # Sort in decreasing dimension
        S_min2 = sorted(simplices, key=lambda x: (len(x), x), reverse=True)
        # Build F_min as a list of pairs (simplex, weight)
        F_min2 = [[s, 0] for s in S_min2]  # F_min: list[list[int], int]
        # Rebuild a dictionary temporarily for sorting
        F_min_dict2 = {tuple(s): w for s, w in F_min2}
        # Re-sort S_min based on F value and dimension
        S_min2 = sorted(S_min2, key=lambda s: (-F_min_dict2[tuple(s)], -len(s)))
        """

        """
        # Group S and F by dimension
        S_by_dim_min = defaultdict(list)
        F_by_dim_min = defaultdict(list)

        for sigma in S_min1:
                S_by_dim_min[len(sigma)-1].append(sigma)
        for sigma, val in F_min1.items():
                F_by_dim_min[len(sigma)-1].append((sigma, val))

        # Prepare numpy arrays for each dimension
        S_arrays_min = {}
        F_arrays_min = {}

        for dim, simplices in S_by_dim_min.items():
                S_arrays_min[dim] = pad_and_pack_S(simplices, dim)

        for dim, simplices_vals in F_by_dim_min.items():
                F_arrays_min[dim] = pad_and_pack_F(simplices_vals, dim)
        
        S_buffer_min = np.vstack([S_arrays_min[dim] for dim in sorted(S_arrays_min)])
        F_buffer_min = np.vstack([F_arrays_min[dim] for dim in sorted(F_arrays_min)])

        # Then pass S_arrays_min[dim] and F_arrays_min[dim] to C++ binding functions
        # for each dimension dim.
        
        for dim in S_arrays_min:
                assert isinstance(S_arrays_min[dim], np.ndarray), f"S_arrays_min[{dim}] is not an ndarray"
                assert isinstance(F_arrays_min[dim], np.ndarray), f"F_arrays_min[{dim}] is not an ndarray"
                assert S_arrays_min[dim].dtype == np.int32, f"S_arrays_min[{dim}] has dtype {S_arrays_min[dim].dtype}"
                assert F_arrays_min[dim].dtype == np.int32, f"F_arrays_min[{dim}] has dtype {F_arrays_min[dim].dtype}"
                assert S_arrays_min[dim].ndim == 2, f"S_arrays_min[{dim}] is not 2D"
                assert F_arrays_min[dim].ndim == 2, f"F_arrays_min[{dim}] is not 2D"
        check_inputs(S_arrays_min, F_arrays_min)
        """

        S_buffer_min = pad_and_pack_S_general(S_min1)
        F_buffer_min = pad_and_pack_F_general(F_min1)
        
        start_min = time.time()
        min, n_crit_min = ms.Min(S_buffer_min, F_buffer_min)
        end_min = time.time()
        time_min = end_min - start_min
        print(f"Min C++/python : {time_min:.6f} secondes")

        print(f"Différence : {time_min2 - time_min:.6f} secondes | Ratio (python pur / C++): {time_min2 / time_min:.2f}\n")
        
#pr = cProfile.Profile()
#pr.enable()

timer_comparison()

"""
pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
"""

# To run the file from the root : python3 -m src.python.python_tests.timer

"""
Pour lancer test.py : 
1) Entrer dans l'environnement virtuel : source venv/bin/activate
2) Lancer test.py : python3 main.py
3) Quitter l'environnement virtuel : deactivate
"""