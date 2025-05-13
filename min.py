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

# computes a minimal increasing Morse sequence obtained 
# from a cosimplicial complex S weighted by a function F.

# Commentaries have been done on max.py
def Min(S, st, F):
    T = dict()
    Sdict = dict()
    U = list()
    MorseSequence=list()
    N = len(S)

    for s in st.simplices():
        T[s] = False
        Sdict[s] = False

    rho = dict()
    for s in S:
        Sdict[s] = True
        nb = nbcoboundary(st, s, Sdict)
        rho[s] = nb
        if nb == 1:
             U.append(s)

    i = 0
    while i<N:
        while U:
            tau = U.pop(0)
            if rho[tau] == 1:
                sigma = next(s for s in coboundary(st, tau, Sdict) if not T[s])
                if F[sigma] == F[tau]:
                    MorseSequence.append([tau, sigma])
                    T[tau] = True
                    T[sigma] = True
                    for mu in boundary(sigma, Sdict)+boundary(tau, Sdict):
                        rho[mu] = rho[mu] - 1
                        if rho[mu] == 1:
                             U.append(mu)
        while i<N and T[S[i]]:
            i += 1
        if i<N:
            sigma = S[i]
            MorseSequence.append([sigma])
            T[sigma] = True
            for tau in boundary(sigma, Sdict):
                    rho[tau] = rho[tau] - 1
                    if rho[tau] == 1:
                        U.append(tau)
    return MorseSequence[::-1]


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
S = sorted(st.simplices(), key=lambda x: (len(x), x))[::-1]
F = dict()
for s in S:
    F[s] = 0
F[(0, 1, 2)] = 0
S = sorted(st.simplices(), key=lambda s: (-F[s], -len(s)))
Min(S, st, F)

#print(f"st = {st.simplices()}")
#print(f"S = {S}")

print(f" Min = {Min(S, st, F)} ")
"""