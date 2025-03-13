from simplextree import SimplexTree

#Reproduis l'union entre un élément elmt et une liste l
def union(elmt, l):
    if elmt not in l:
        l.append(elmt)
        

#Retourne la dimension d'un simplexe
def dim(simplexe):
    return len(simplexe) - 1


#Retourne le bord d'un simplexe
def bord(simplexe):
    n = dim(simplexe)
    st = SimplexTree([simplexe])
    bord = [sigma for sigma in st.simplices(n-1)]
    return bord


#Retourne le cobord d'un simplexe
def cobord(simplexe, K_init):
    cobord = K_init.cofaces(simplexe)
    i = 0  # Index pour parcourir la liste cobord
    while i < len(cobord):  # Equivalent à "do...while" avec une seule condition
        elmt = cobord[i]
        if dim(elmt) != dim(simplexe) + 1:
            cobord.remove(elmt)
        else:
            i += 1  # Passer à l'élément suivant uniquement si aucun élément n'est supprimé
    return cobord


#Retourne le simplexe v tel que :
#   - v est dans bord
#   - v n'est pas dans S
def find_out(bord, S): 
    possibilities = []
    #print(f"bord = {bord}")
    #print(f"S = {S}")
    for v in bord:
        #print(f"elmt de bord : {v}")
        if v not in S:
            possibilities.append(v)
    return possibilities[0]

#Trie la liste L de simplexes par dimension décroissante
def tri_dim(L):
    L.sort(key=dim, reverse=True)
    
# =============================================================================================================================================================================== #

def morse_seq_croissante(K_init):
    
    # Trier les simplexes de K_init par dimensions croissante
    K = K_init.simplices()
    n = len(K)
    
    # Initialisation des variables
    i = 0
    S = list()
    W = list()
    L = list()
    ro = {s: 0 for s in K}
    n_crit = 0
    
    while i < n-1:  # i va de 0 à n-1
        sigma = K[i]
        S.append(sigma)
        W.append(sigma) 
        n_crit += 1
        
        for tau in cobord(sigma, K_init): 
            ro[tau] += 1
            L.append(tau)
            tri_dim(L)
        
        while len(L) > 0:
            tau = L.pop()
            
            if ro[tau] == dim(tau):
                v = find_out(bord(tau), S) 
                W.append((v, tau))
                union(v, S)
                union(tau, S)
                
                c1 = cobord(v, K_init)
                c2 = cobord(tau, K_init)
                cobords = []
                for elmt in c1 + c2:
                    union(elmt, cobords)
                for mu in cobords: 
                    ro[mu] += 1
                    union(mu, L)
                tri_dim(L)
                    
        while K[i] in S and i < n-1: 
            i += 1
            
    return W, n_crit
    
# =============================================================================================================================================================================== #

#K_init = SimplexTree([[1,2,3]])

K_init = SimplexTree([[1, 5, 7], [1, 2, 7],    # Haut gauche
                        [2, 7, 9], [2, 3, 9],  # Haut milieu
                        [3, 5, 9], [1, 3, 5],  # Haut droit
                        [5, 4, 6], [5, 6, 7],  # Milieu gauche
                        [7, 6, 8], [7, 8, 9],  # Milieu centre
                        [9, 8, 4], [9, 4, 5],  # Milieu droit
                        [1, 2, 4], [2, 4, 6],  # Bas gauche
                        [2, 3, 6], [3, 6, 8],  # Bas milieu
                        [1, 3, 8], [1, 4, 8]])  # Bas droit

seq, n_crit = morse_seq_croissante(K_init)
