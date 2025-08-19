from simplextree import SimplexTree
import morse_sequence._core
from morse_sequence._core import MorseSequence

st = SimplexTree([[1,2,3]])
cpp_ms = morse_sequence._core.cpp_ms_from_py_st(st)

print("python st type", type(st))
print("cpp st type", type(cpp_ms))

ms = MorseSequence(st)
#print(ms.simplices())
W, n_crit = ms.increasing()


"""
from typing import Union, Tuple, List

def show_morse_frame(title: str, mf: dict, W: list) -> None:
    print(f"{title} ({len(mf)} entrées):\n")

    for item in W:
        if isinstance(item, (list, tuple)) and all(isinstance(x, int) for x in item):
            key = tuple(sorted(item))
            val = mf.get(key, [])
            print(f"Key (Critical simplex): {key} -> Value:", end=" ")
            for v in val:
                print(v if v is not None else "None", end=" ")
            print()

        elif isinstance(item, (list, tuple)) and len(item) == 2:
            lower = tuple(sorted(item[0]))
            upper = tuple(sorted(item[1]))
            val_lower = mf.get(lower, [])
            val_upper = mf.get(upper, [])

            print("Pair of simplices:")
            print(f"  Key (Lower pair) {lower} -> Value:", end=" ")
            for v in val_lower:
                print(v if v is not None else "None", end=" ")
            print()

            print(f"  Key (Upper pair) {upper} -> Value:", end=" ")
            for v in val_upper:
                print(v if v is not None else "None", end=" ")
            print()
        print()

    print()



st = SimplexTree([[1, 5, 7], [1, 2, 7],    # Top left
                 [2, 7, 9], [2, 3, 9],    # Top middle
                 [3, 5, 9], [1, 3, 5],    # Top right
                 [5, 4, 6], [5, 6, 7],    # Middle left
                 [7, 6, 8], [7, 8, 9],    # Middle center
                 [9, 8, 4], [9, 4, 5],    # Middle right
                 [1, 2, 4], [2, 4, 6],    # Bottom left
                 [2, 3, 6], [3, 6, 8],    # Bottom middle
                 [1, 3, 8], [1, 4, 8]])   # Bottom right

ms = MorseSequence(st)
W, n_crit = ms.increasing()

reference_map = ms.reference_map(W)
coreference_map = ms.coreference_map(W)

print("W (Morse sequence) :\n", W, "\n")
show_morse_frame("Reference map",   reference_map, W)
show_morse_frame("Coreference map", coreference_map, W)
"""
# To run the file from the root : python3 -m tests.python_tests.test_reference