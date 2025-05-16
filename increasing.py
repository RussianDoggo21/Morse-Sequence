from simplextree import SimplexTree
from itertools import combinations

"""
# Return the boundary of a simplex
def boundary(simplex):
    n = dim(simplex)
    st = SimplexTree([simplex])
    boundary = [sigma for sigma in st.simplices(n-1)]
    return boundary
"""

"""
# Return the coboundary of a simplex
def coboundary(simplexe, K_init):
    coboundary = K_init.cofaces(simplexe)
    i = 0  # Index to browse the list coboundary 
    while i < len(coboundary):  # Equivalent to "do...while" with an only condition
        elmt = coboundary[i]
        if dim(elmt) != dim(simplexe) + 1:
            coboundary.remove(elmt)
        else:
            i += 1  # Go to the next element only if no element has been removed
    return coboundary
"""

# Return the dimension of a simplex
def dim(simplex):
    return len(simplex) - 1

# Compute the boundary of the simplexe sigma in the complex S
def boundary(sigma):
    if len(sigma) > 1:
        return [tuple(s) for s in combinations(sigma,len(sigma)-1)]
    return list()

# Compute the coboundary of the simplexe sigma in the complex S
def coboundary(sigma, st): 
    return [s for s in st.cofaces(sigma) if (len(s) == len(sigma) + 1)]

# Compute the length of the coboundary of the simplexe sigma (its number of cofaces) in the complex S
def nbcoboundary(sigma, st): 
    return len(coboundary(sigma, st))

# Compute the length of the boundary of the simplexe sigma (its number of faces) in the complex S
def nbboundary(st, sigma): 
    return len(boundary(sigma))

# Return the simplex v such that :
#   - v is in s_list
#   - v is not in S
def find_out(s_list, S): 
    possibilities = []
    for v in s_list:
        if v not in S:
            possibilities.append(v)
    return possibilities[0]

# Sort the simplices list L by decreasing dimension 
def tri_dim(L):
    L.sort(key=dim, reverse=True)
    
# Main function
def morse_seq_increasing(K_init):
    
    # Sort the simplices of K_init by increasing dimension 
    K = K_init.simplices()
    n = len(K)
    
    # Initialization of the variables
    i = 0 # Index allowing to browse the simplices of K
    W = list() # Morse Sequence to be returned
    S = list() # Contains the already "used" simplices in the construction of W
    L = list() # Contains the simplices of dimension p that can form a free pair (p-1, p) for W
    rho = {s: 0 for s in K} # Contains the number of faces already "used" for W of each simplex s in K
    n_crit = 0 # Counter of critical simplices in W
    
    while i < n-1:  # i goes from 0 to n-1
        
        # Here, sigma is a critical simplex
        sigma = K[i]
        S.append(sigma)
        W.append(sigma) 
        n_crit += 1
        
        # Update of rho and L
        for tau in coboundary(sigma, K_init): 
            rho[tau] += 1
            L.append(tau)
            tri_dim(L) # Sort L by decreasing dimension
        
        # While we can add free pairs to W
        while len(L) > 0:
            tau = L.pop() # Possible first element of the free pair, dimension p
            if rho[tau] == dim(tau): # Check if tau can be used for a free pair of W
                v = find_out(boundary(tau), S) # Second element of the free pair, dimension p-1
                
                # Update of W and S
                W.append((v, tau))
                S.append(v)
                S.append(tau)
                
                # Update of rho and L
                c1 = coboundary(v, K_init) # all simplices of dimension p
                c2 = coboundary(tau, K_init ) # all simplices of dimension p+1
                coboundarys = []
                for elmt in c1 + c2:
                    coboundarys.append(elmt) # No risk of duplicate due to the difference in dimension
                for mu in coboundarys: 
                    rho[mu] += 1
                    L.append(mu)
                tri_dim(L) # Sort L by decreasing dimension
        
        # If the current simplex has already been "used" for W and we didn't
        while K[i] in S and i < n-1: 
            i += 1
            
    return W, n_crit
    
# =============================================================================================================================================================================== #
import numpy as np
#K_init = SimplexTree([[1,2,3]])
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

K_init = SimplexTree([[1, 5, 7], [1, 2, 7],    # Haut gauche
                        [2, 7, 9], [2, 3, 9],  # Haut milieu
                        [3, 5, 9], [1, 3, 5],  # Haut droit
                        [5, 4, 6], [5, 6, 7],  # Milieu gauche
                        [7, 6, 8], [7, 8, 9],  # Milieu centre
                        [9, 8, 4], [9, 4, 5],  # Milieu droit
                        [1, 2, 4], [2, 4, 6],  # Bas gauche
                        [2, 3, 6], [3, 6, 8],  # Bas milieu
                        [1, 3, 8], [1, 4, 8]])  # Bas droit

seq, n_crit = morse_seq_increasing(K_init)

print(seq, "\n")
print(n_crit)
