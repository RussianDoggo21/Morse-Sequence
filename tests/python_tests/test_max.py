from simplextree import SimplexTree
from src.python.max import Max

st = SimplexTree([[1, 5, 7], [1, 2, 7],    # Top left
                 [2, 7, 9], [2, 3, 9],    # Top middle
                 [3, 5, 9], [1, 3, 5],    # Top right
                 [5, 4, 6], [5, 6, 7],    # Middle left
                 [7, 6, 8], [7, 8, 9],    # Middle center
                 [9, 8, 4], [9, 4, 5],    # Middle right
                 [1, 2, 4], [2, 4, 6],    # Bottom left
                 [2, 3, 6], [3, 6, 8],    # Bottom middle
                 [1, 3, 8], [1, 4, 8]])   # Bottom right


S_max = sorted(st.simplices(), key=lambda x: (len(x), x))
F_max = dict()
for s in S_max:
    F_max[s] = 0
F_max[(0, 1, 2)] = 1
S_max = sorted(st.simplices(), key=lambda s: (F_max[s], len(s)))

"""
# Convert all simplices to list[int] instead of tuple (used in Max and Min)
simplices = [list(s) for s in st.simplices()]

# Sort in increasing dimension
S_max = sorted(simplices, key=lambda x: (len(x), x))

# Build F_max as a list of pairs (simplex, weight)
F_max = [[s, 0] for s in S_max]  # F_max: list[list[int], int]

# Rebuild a dictionary temporarily for sorting
F_max_dict = {tuple(s): w for s, w in F_max}

# Re-sort S_max based on F value and dimension
S_max = sorted(S_max, key=lambda s: (F_max_dict[tuple(s)], len(s)))
"""
max, n_crit = Max(S_max, st, F_max)
print(f" n_crit = {n_crit}\n Max = {max} ")

# To run the file from the root : python3 -m tests.python_tests.test_max