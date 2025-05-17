from simplextree import SimplexTree
from itertools import combinations

# Compute the boundary of the simplex sigma in the complex S
def boundary(sigma, S):
    if len(sigma) > 1:
        return [tuple(s) for s in combinations(sigma, len(sigma) - 1) if S[s]]
    return list()

# Compute the coboundary of the simplex sigma in the complex S
def coboundary(st, sigma, S): 
    return [s for s in st.cofaces(sigma) if (len(s) == len(sigma) + 1) and S[s]]

# Compute the number of cofaces of the simplex sigma in the complex S
def nbcoboundary(st, sigma, S): 
    return len(coboundary(st, sigma, S))

# Compute the number of faces of the simplex sigma in the complex S
def nbboundary(st, sigma, S): 
    return len(boundary(sigma, S))


# =============================================================================================================================================================================== #

# Computes a minimal increasing Morse sequence obtained 
# from a cosimplicial complex S weighted by a function F.
def Min(S, st, F):

    T = dict()          # Boolean dictionary: T[s] is True if simplex s has been used in the Morse Sequence
    Sdict = dict()      # Boolean dictionary: Sdict[s] is True if simplex s belongs to the input cosimplicial complex S
    U = list()          # List of simplices sigma such that sigma has only one coface tau: (sigma, tau) is a free pair
    MorseSequence = []  # Morse Sequence to be returned
    N = len(S)          # Number of simplices in S
    rho = dict()        # Dictionary: rho[s] is the number of cofaces of s
    n_crit = 0          # Counter of critical simplices in the Morse Sequence
    i = 0  # Index to browse the simplices of S

    # Initialization of dictionaries T and Sdict
    for s in st.simplices():
        T[s] = False       # No simplex is used at the beginning
        Sdict[s] = False   # Initialize Sdict to False for all simplices in the full complex

    # Initialization of rho and U
    for s in S:
        Sdict[s] = True
        nb = nbcoboundary(st, s, Sdict)  # Count the number of cofaces of s in S
        rho[s] = nb
        if nb == 1:
            U.append(s)  # Candidate for a free pair (sigma, tau)

    # While we haven't used all simplices in S 
    while i < N:

        # While we can still form free pairs (sigma, tau)
        while U:
            sigma = U.pop(0)
            if rho[sigma] == 1:
                # There should be only one such tau with T[tau] == False
                tau = next(s for s in coboundary(st, sigma, Sdict) if not T[s])
                
                # Check if (tau, sigma) can form a free pair according to F
                if F[sigma] == F[tau]:

                    # Update of the Morse Sequence and usage status
                    MorseSequence.append([sigma, tau])
                    T[tau] = True
                    T[sigma] = True

                    # Update of rho and U
                    for mu in boundary(sigma, Sdict) + boundary(tau, Sdict):
                        rho[mu] -= 1
                        if rho[mu] == 1:
                            U.append(mu)

        # Skip all simplices already used in the Morse Sequence
        while i < N and T[S[i]]:
            i += 1

        # If we still have unused simplices, add a critical one
        if i < N:
            sigma = S[i]
            MorseSequence.append([sigma])  # Add as a critical simplex
            n_crit += 1
            T[sigma] = True

            # Update of rho and U
            for tau in boundary(sigma, Sdict):
                rho[tau] -= 1
                if rho[tau] == 1:
                    U.append(tau)

    return MorseSequence[::-1], n_crit  # Reverse the list to obtain an increasing sequence



# =============================================================================================================================================================================== #

st = SimplexTree([[1, 5, 7], [1, 2, 7],    # Haut gauche
                        [2, 7, 9], [2, 3, 9],  # Haut milieu
                        [3, 5, 9], [1, 3, 5],  # Haut droit
                        [5, 4, 6], [5, 6, 7],  # Milieu gauche
                        [7, 6, 8], [7, 8, 9],  # Milieu centre
                        [9, 8, 4], [9, 4, 5],  # Milieu droit
                        [1, 2, 4], [2, 4, 6],  # Bas gauche
                        [2, 3, 6], [3, 6, 8],  # Bas milieu
                        [1, 3, 8], [1, 4, 8]])  # Bas droit
S = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
F = dict()
for s in S:
    F[s] = 0
F[(0, 1, 2)] = 0
S = sorted(st.simplices(), key=lambda s: (-F[s], -len(s)))
Min(S, st, F)

#print(f"st = {st.simplices()}")
#print(f"S = {S}")
min, n_crit = Min(S, st, F)
print(f"n_crit = {n_crit}\n Min = {min} ")
