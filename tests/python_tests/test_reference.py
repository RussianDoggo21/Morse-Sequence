from morse_sequence import MorseSequence
from simplextree import SimplexTree

def show_morse_frame(title: str, mf: dict) -> None:
    """Show a Python morse_frame {Key(tuple) : Value [list]}."""
    print(f"{title} ({len(mf)} Entries) :\n")
    for idx, (key, val) in enumerate(mf.items(), 1):
        print(f"{idx:2d}. key = {key}  →  value = {val}")
    print()                 # ligne vide finale

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
W, n_crit = ms.increasing(st)

reference_map = ms.reference_map(W)
coreference_map = ms.coreference_map(W)

print("W (Morse sequence) :\n", W, "\n")
show_morse_frame("Reference map",   reference_map)
show_morse_frame("Coreference map", coreference_map)

# To run the file from the root : python3 -m morse_sequence.tests.test_reference