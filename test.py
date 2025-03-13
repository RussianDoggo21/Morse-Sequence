from simplextree import SimplexTree
from morse_sequence import MorseSequence

# Créer une instance de SimplexTree
tree = SimplexTree([[1, 5, 7], [1, 2, 7], [2, 7, 9], [2, 3, 9], [3, 5, 9], [1, 3, 5], [5, 4, 6], [5, 6, 7], [7, 6, 8],
                   [7, 8, 9], [9, 8, 4], [9, 4, 5], [1, 2, 4], [2, 4, 6], [2, 3, 6], [3, 6, 8], [1, 3, 8], [1, 4, 8]])


# Créer une instance de MorseSequence
ms = MorseSequence(tree)

# Appeler les méthodes
ms_dec, n_crit = ms.morse_seq_dec(tree)
ms_crois, n_crit2 = ms.morse_seq_crois(tree)


print(f"ms_dec :\n {ms_dec},\n n_crit = {n_crit}\n\n")
print(f"ms_crois :\n {ms_crois},\n n_crit = {n_crit2}\n\n")
