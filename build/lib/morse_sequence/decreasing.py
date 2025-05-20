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
def coboundary(simplex, K_init):
    coboundary = K_init.cofaces(simplex)
    i = 0  # Index to browse the list coboundary
    while i < len(coboundary):  # Equivalent to "do...while" with only one condition
        elmt = coboundary[i]
        if dim(elmt) != dim(simplex) + 1:
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
def nbcoboundary(st, sigma, S): 
    return len(coboundary(st, sigma, S))

# Compute the length of the boundary of the simplexe sigma (its number of faces) in the complex S
def nbboundary(st, sigma, S): 
    return len(boundary(sigma, S))

# Return the simplex v such that:
#   - v is in s_list
#   - v is not in S
def find_out(s_list, S): 
    possibilities = []
    for v in s_list:
        if v not in S:
            possibilities.append(v)
    if len(possibilities) > 0:
        return possibilities[0]
    return None


# Main function for decreasing Morse sequence
def morse_seq_decreasing(K_init):

    # Get the simplices of K_init sorted by decreasing dimension
    K = K_init.simplices()
    K.sort(key=lambda s: dim(s), reverse=True)
    n = len(K)
    
    # Initialization of variables
    i = 0 # Index to browse the simplices of K
    W = list() # Morse sequence to be returned
    S = list() # Contains already "used" simplices for W
    rho = {s: len(coboundary(s, K_init)) for s in K} # Number of cofaces currently available of each simplex
    L = [s for s in K if rho[s] == 1] # Simplices of dimension p that can potentially form a free pair (p, p+1) for W
    n_crit = 0 # Counter of critical simplices
    
    # Main loop
    while i < n:
        # While we can add free pairs to W
        while len(L) > 0:
            v = L.pop() # Possible first element of the free pair, p dimension
            tau = find_out(coboundary(v, K_init), S) # Possible second element of the free pair, p+1 dimension
            if tau is None: # v = None means that (v, tau) isn't a free pair for W
                continue # We go to the next v in L

            W.append((v, tau)) # We add (v,tau) in W if it is  a free pair
            
            # Update S
            S.append(v)
            S.append(tau)
            
            # Update rho and L
            c1 = boundary(v)
            c2 = boundary(tau)
            boundaries = []
            for elmt in c1 + c2:
                if elmt not in boundaries:
                    boundaries.append(elmt)
            for mu in boundaries:
                rho[mu] -= 1
                if rho[mu] == 1 and mu not in L:
                    L.append(mu)
        
        # Skip already used simplices
        while i < n and K[i] in S:
            i += 1
        if i == n:
            return W, n_crit
        
        # If we reach a critical simplex
        critical = K[i]
        W.append([critical])
        n_crit += 1
        
        S.append(critical)
        
        # Update rho and L with the boundary of the critical simplex
        for simplex in boundary(critical):
            rho[simplex] -= 1
            if rho[simplex] == 1 and simplex not in L:
                L.append(simplex)
    
    return W, n_crit

            



