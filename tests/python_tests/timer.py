# File to measure the time complexity of each of the 4 algorithms
# In both pure python and hybrid code

from simplextree import SimplexTree
from morse_sequence import MorseSequence

from src.python.max import Max
from src.python.min import Min
from src.python.increasing import morse_seq_increasing
from src.python.decreasing import morse_seq_decreasing

import numpy as np
import time
import csv
import os

CSV_FILE = "/home/kali/Téléchargements/Morse-Sequence/tests/graphs/timing_results_python.csv"

def write_timing_result(size, operation, language, time_seconds):
    header = ["size", "operation", "language", "time_seconds"]
    write_header = not os.path.exists(CSV_FILE)
    with open(CSV_FILE, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(header)
        writer.writerow([size, operation, language, round(time_seconds, 6)])


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

def timer_comparison():
    total_start = time.time()
    list_faces = [10, 20, 50, 60, 75, 100]
    for k in list_faces:
        print(f"\n======================= Grid case {k} x {k} =======================\n")

        # --- Block: Create SimplexTree + MorseSequence ---
        start_init = time.time()
        st = SimplexTree(MakeFacesVectorized1(k, k))
        ms = MorseSequence(st)
        end_init = time.time()
        print(f"Initialization (SimplexTree + MorseSequence): {end_init - start_init:.6f} seconds\n")
        write_timing_result(k, "init", "python", end_init - start_init)

        # size,operation,language,time_seconds

        # ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        start_dec2 = time.time()
        ms_dec2, n_crit_dec2 = morse_seq_decreasing(st)
        end_dec2 = time.time()
        time_dec2 = end_dec2 - start_dec2
        print(f"Decreasing pure Python: {time_dec2:.6f} seconds")
        write_timing_result(k, "decreasing", "python", time_dec2)

        start_dec = time.time()
        ms_dec, n_crit_dec = ms.decreasing(st)
        end_dec = time.time()
        time_dec = end_dec - start_dec
        print(f"Decreasing C++/Python: {time_dec:.6f} seconds")
        write_timing_result(k, "decreasing", "cpp/python", time_dec)

        print(f"Difference: {time_dec2 - time_dec:.6f} seconds | Ratio (pure Python / C++): {time_dec2 / time_dec:.2f}\n")
        
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        start_crois2 = time.time()
        ms_crois2, n_crit_crois2 = morse_seq_increasing(st)
        end_crois2 = time.time()
        time_crois2 = end_crois2 - start_crois2
        print(f"Increasing pure Python: {time_crois2:.6f} seconds")
        write_timing_result(k, "increasing", "python", time_crois2)

        start_crois = time.time()
        ms_crois, n_crit_crois = ms.increasing(st)
        end_crois = time.time()
        time_crois = end_crois - start_crois
        print(f"Increasing C++/Python: {time_crois:.6f} seconds")
        write_timing_result(k, "increasing", "cpp/python", time_crois)

        print(f"Difference: {time_crois2 - time_crois:.6f} seconds | Ratio (pure Python / C++): {time_crois2 / time_crois:.2f}\n")
        
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        start_max_pre = time.time()
        S_max1 = sorted(st.simplices(), key=lambda x: (len(x), x))
        F_max1 = {s: 0 for s in S_max1}
        S_max1 = sorted(st.simplices(), key=lambda s: (F_max1[s], len(s)))
        end_max_pre = time.time()
        print(f"Max preprocessing pure Python: {end_max_pre - start_max_pre:.6f} seconds")
        write_timing_result(k, "max_preprocessing", "python", end_max_pre - start_max_pre)

        start_max2 = time.time()
        max2 = Max(S_max1, st, F_max1)
        end_max2 = time.time()
        time_max2 = end_max2 - start_max2
        print(f"Max pure Python: {time_max2:.6f} seconds")
        write_timing_result(k, "max", "python", time_max2)

        start_pack_max = time.time()
        S_buffer_max = pad_and_pack_S_general(S_max1)
        F_buffer_max = pad_and_pack_F_general(F_max1)
        end_pack_max = time.time()
        print(f"Max packing for C++: {end_pack_max - start_pack_max:.6f} seconds")
        write_timing_result(k, "max_preprocessing", "cpp/python", end_pack_max - start_pack_max)

        start_max = time.time()
        max, n_crit_max = ms.Max(S_buffer_max, F_buffer_max)
        end_max = time.time()
        time_max = end_max - start_max
        print(f"Max C++/Python: {time_max:.6f} seconds")
        write_timing_result(k, "max", "cpp/python", time_max)

        print(f"Difference: {time_max2 - time_max:.6f} seconds | Ratio (pure Python / C++): {time_max2 / time_max:.2f}\n")

        # ------------------------------------------------------------------------------------------------------------------------------------------------------------ #
        start_min_pre = time.time()
        S_min1 = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
        F_min1 = {s: 0 for s in S_min1}
        S_min1 = sorted(st.simplices(), key=lambda s: (-F_min1[s], -len(s)))
        end_min_pre = time.time()
        print(f"Min preprocessing pure Python: {end_min_pre - start_min_pre:.6f} seconds")
        write_timing_result(k, "min_preprocessing", "python", end_min_pre - start_min_pre)

        start_min2 = time.time()
        min2 = Min(S_min1, st, F_min1)
        end_min2 = time.time()
        time_min2 = end_min2 - start_min2
        print(f"Min pure Python: {time_min2:.6f} seconds")
        write_timing_result(k, "min", "python", time_min2)

        start_pack_min = time.time()
        S_buffer_min = pad_and_pack_S_general(S_min1)
        F_buffer_min = pad_and_pack_F_general(F_min1)
        end_pack_min = time.time()
        print(f"Min packing for C++: {end_pack_min - start_pack_min:.6f} seconds")
        write_timing_result(k, "min_preprocessing", "cpp/python", end_pack_min - start_pack_min)

        start_min = time.time()
        min, n_crit_min = ms.Min(S_buffer_min, F_buffer_min)
        end_min = time.time()
        time_min = end_min - start_min
        print(f"Min C++/Python: {time_min:.6f} seconds")
        write_timing_result(k, "min", "cpp/python", time_min)

        print(f"Difference: {time_min2 - time_min:.6f} seconds | Ratio (pure Python / C++): {time_min2 / time_min:.2f}\n")

    total_end = time.time()
    print(f"\n⏱️ Total time timer_comparison(): {total_end - total_start:.2f} seconds")
    write_timing_result(0, "total", "python", total_end - total_start)




timer_comparison()

# To run the file from the root : python3 -m tests.python_tests.timer

"""
Pour lancer test.py : 
1) Entrer dans l'environnement virtuel : source venv/bin/activate
2) Lancer test.py : python3 main.py
3) Quitter l'environnement virtuel : deactivate
"""