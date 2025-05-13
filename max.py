from simplextree import SimplexTree
from itertools import combinations

# Compute the boundary of the simplexe sigma in the complex S
def boundary(sigma, S):
    if len(sigma) > 1:
        return [tuple(s) for s in combinations(sigma,len(sigma)-1) if S[s]]
    return list()

# Compute the coboundary of the simplexe sigma in the complex S
def coboundary(st, sigma, S): 
    return [s for s in st.cofaces(sigma) if (len(s) == len(sigma) + 1) and S[s]]

# Compute the length of the coboundary of the simplexe sigma (its number of cofaces) in the complex S
def nbcoboundary(st, sigma, S): 
    return len(coboundary(st, sigma, S))

# Compute the length of the boundary of the simplexe sigma (its number of faces) in the complex S
def nbboundary(st, sigma, S): 
    return len(boundary(sigma, S))

def difflist(list1, list2): 
    # Can be used to obtain the differences of two simplicial complexes
    ## st1 = SimplexTree([[0,1,2]])
    ## st2 = SimplexTree([[1,2]])
    ## S = [[0], [0,1], [0,2], [0,1,2]]
    ## assert S == difflist(st1.simplices(), st2.simplices())

    s1 = {tuple(i) for i in list1}
    s2 = {tuple(i) for i in list2}
    l = list(list(i) for i in (s1 - s2))
    return sorted(l, key=len)

# =============================================================================================================================================================================== #

# computes a maximal increasing Morse sequence obtained 
# from a cosimplicial complex S weighted by a function F.

def Max(S, st, F):
    T = dict() # Boolean dictionnary : if T[s] == False, s is still "available to use" for the Morse Sequence
    Sdict = dict() # ??????
    U = list() # List containing simplices v with only one coface tau : (v, tau) is a free pair for the Morse Sequence
    MorseSequence=list() # Morse Sequence to return 
    N = len(S) # Number of simplices in S (and st)

    for s in st.simplices():
        T[s] = False # No simplices has been "used" yet
        Sdict[s] = False # ??

    rho = dict() # Dictionnary : if rho[s] == n, it means s has n faces
    for s in S:
        Sdict[s] = True # ???
        nb = nbboundary(st, s, Sdict)
        rho[s] = nb
        if nb == 1:
             U.append(s)

    i = 0
    while i<N: # While we haven't used all simplices in S (and st)

        while U: # While we can still add free pairs
            tau = U.pop(0) 
            if rho[tau] == 1: # First verification on tau
                sigma = next(s for s in boundary(tau, Sdict) if not T[s]) # boundary(tau, Sdict) should return only a single simplex
                                                                          # The condition "if not T[s]" allows us to work with a changing simplicial complex
                if F[sigma] == F[tau]: # Second verification on tau and sigma : is it really a free pair (by the definition of F)
                    
                    # Update of MorseSequence and T
                    MorseSequence.append([sigma, tau]) 
                    T[tau] = True # tau has been "used" in the Morse Sequence
                    T[sigma] = True # sigma has been "used" in the Morse Sequence
                    
                    # Update of rho and then of U
                    for mu in coboundary(st, sigma, Sdict)+coboundary(st, tau, Sdict):
                        rho[mu] = rho[mu] - 1
                        if rho[mu] == 1:
                             U.append(mu)

        # We skip all the simplices that have already been "used" in the Morse Sequence
        while i<N and T[S[i]]:
            i += 1
        
        # Now that we have added all free pairs (loop while U), the next simplice should be a critical one
        # Update of MorseSequence, T, rho and U accordingly
        if i<N:
            sigma = S[i]
            MorseSequence.append([sigma])
            T[sigma] = True
            for tau in coboundary(st, sigma, Sdict):
                    rho[tau] = rho[tau] - 1
                    if rho[tau] == 1:
                        U.append(tau)

    return MorseSequence

# =============================================================================================================================================================================== #
"""
st = SimplexTree([[1, 5, 7], [1, 2, 7],    # Haut gauche
                        [2, 7, 9], [2, 3, 9],  # Haut milieu
                        [3, 5, 9], [1, 3, 5],  # Haut droit
                        [5, 4, 6], [5, 6, 7],  # Milieu gauche
                        [7, 6, 8], [7, 8, 9],  # Milieu centre
                        [9, 8, 4], [9, 4, 5],  # Milieu droit
                        [1, 2, 4], [2, 4, 6],  # Bas gauche
                        [2, 3, 6], [3, 6, 8],  # Bas milieu
                        [1, 3, 8], [1, 4, 8]])  # Bas droit
S = sorted(st.simplices(), key=lambda x: (len(x), x))
F = dict()
for s in S:
    F[s] = 0
F[(0, 1, 2)] = 1
S = sorted(st.simplices(), key=lambda s: (F[s], len(s)))
#print(f"st = {st.simplices()}")
#print(f"S = {S}")
print(f" Max = {Max(S, st, F)} ")
"""