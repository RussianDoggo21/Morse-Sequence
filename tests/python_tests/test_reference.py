from simplextree import SimplexTree
from morse_sequence._core import MorseSequence

def show_morse_frame(title: str, mf: list) -> None:
    print(f"{title} ({len(mf)} entries):\n")
    for entry in mf:
        if not isinstance(entry, tuple) or len(entry) != 2:
            print(f"Invalid entry, skipping: {entry}")
            continue
        key, vals = entry
        key_set = set(key)
        if len(vals) == 0:
            print(f"Key (Critical simplex): {key_set} -> Value: {{}}")
        else:
            print(f"Key: {key_set} -> Values: {', '.join(str(set(v)) for v in vals)}")
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
#coreference_map = ms.coreference_map(W)

#print("W (Morse sequence) :\n", W, "\n")
print("Raw reference_map:", reference_map)
#print("Raw co-reference_map:", coreference_map)

#show_morse_frame("Reference map",   reference_map)
#show_morse_frame("Coreference map", coreference_map)


# To run the file from the root : python3 -m tests.python_tests.test_reference