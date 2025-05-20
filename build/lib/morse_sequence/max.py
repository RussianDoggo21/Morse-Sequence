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

# =============================================================================================================================================================================== #

# computes a maximal increasing Morse sequence obtained 
# from a cosimplicial complex S weighted by a function F.

def Max(S, st, F):

    T = dict() # Boolean dictionnary : if T[s] == False, s is still "available to use" for the Morse Sequence
    Sdict = dict() # Boolean dictionnary : if Sdict[s] == True, s is in S
    U = list() # List containing simplices tau with only one face v : (v, tau) is a free pair for the Morse Sequence
    MorseSequence=list() # Morse Sequence to return 
    N = len(S) # Number of simplices in S 
    rho = dict() # Dictionnary : if rho[s] == n, it means s has n faces
    i = 0 # Indew to browse the simplices in S
    n_crit = 0 # Counts the number of critical simplices in W

    for s in st.simplices():
        T[s] = False # No simplices has been "used" yet
        Sdict[s] = False # At first approach, s being in st doesn't mean that s is in S

    for s in S:
        Sdict[s] = True # s is indeed in S
        nb = nbboundary(st, s, Sdict)
        rho[s] = nb
        # if nb == 1, s verifies the condition to be in U
        if nb == 1:
             U.append(s) 

    while i<N: # While we haven't used all simplices in S 

        while U: # While we can still add free pairs
            tau = U.pop(0) # Possible first element of the pair
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
            n_crit += 1
            T[sigma] = True
            for tau in coboundary(st, sigma, Sdict):
                    rho[tau] = rho[tau] - 1
                    if rho[tau] == 1:
                        U.append(tau)

    return MorseSequence, n_crit


