from simplextree import SimplexTree

#Reproduis l'union entre un élément elmt et une liste l

#PRBLM : TROUVE elmt DANS l ALORS QUE CE N'EST PAS CENSE ETRE LE CAS 
#PROBLEME DES CROCHETS exemple : différence entre (2,3) et [(2,3)]
def union(elmt, l):
    if elmt not in l:
        l.append(elmt)
    else :
        print(f"{elmt} est dans la liste") #Plus besoin de union ?
        

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
#   - v est dans cobord
#   - v n'est pas dans S
def find_out2(cobord, S): 
    possibilities = []
    for v in cobord:
        if v not in S:
            possibilities.append(v)
    """ 
    print(f"cobord : {cobord}")
    print(f"S : {S}")
    print(f"possibilities : {possibilities}")
    """
    if len(possibilities) > 0:
        return possibilities[0]
    return None


# ================================================================================================================================================== #

# AVEC LES BONNES MODIFS (+ print à faire)

def morse_seq_decroissante(K_init):
    # Récupérer les simplexes de K_init 
    K = K_init.simplices()
    K.sort(key=lambda s: dim(s), reverse=True)

    n = len(K)
    
    # Initialisation des variables
    i = 0
    S = list()
    W = list()
    ro = {s: len(cobord(s, K_init)) for s in K}
    L = [s for s in K if ro[s]==1]
    n_crit = 0
    crit = []
    
    while i < n: #modif 5, utilisation de i à la place de len(S)
        while len(L) > 0:
            tau = L.pop()
            #print(f"tau = {tau}")
            v = find_out2(cobord(tau, K_init), S)
            if v == None:
                #print(f"No free pair for {tau}\n")
                continue
            #print(f"v = {v}")
            W.append((v, tau))
            #print(f"W = {W}")
            union(v, S)
            union(tau, S)
            c1 = bord(v)
            c2 = bord(tau)
            bords = []
            for elmt in c1 + c2:
                union(elmt, bords)
            for mu in bords: 
                ro[mu] -= 1
                if ro[mu] == 1:
                    union(mu, L)
            
        while i < n and K[i] in S:  #modif 1 (inversion des conditions)
            i+=1
        if i == n:   
            return W, n_crit #modif 2 (ajout du return)
        
        critical = K[i]
        #print(f"Critical = {critical}")
        W.append([critical]) #modif 3 (ajout des '[]')
        n_crit += 1
        #print(f"W = {W}")
        union(critical, S)
        
        for simplex in bord(critical):
            ro[simplex] -= 1
            if ro[simplex] == 1:
                union(simplex, L)
                
        #if len(S) == n:
        if i == n-1:
            print(f"len(S) == n\n")
            print("=================== End of the programm ===================")
            
# ================================================================================================================================================== #

K_init = SimplexTree([[1, 5, 7], [1, 2, 7],    # Haut gauche
                        [2, 7, 9], [2, 3, 9],  # Haut milieu
                        [3, 5, 9], [1, 3, 5],  # Haut droit
                        [5, 4, 6], [5, 6, 7],  # Milieu gauche
                        [7, 6, 8], [7, 8, 9],  # Milieu centre
                        [9, 8, 4], [9, 4, 5],  # Milieu droit
                        [1, 2, 4], [2, 4, 6],  # Bas gauche
                        [2, 3, 6], [3, 6, 8],  # Bas milieu
                        [1, 3, 8], [1, 4, 8]])  # Bas droit

#K_init = SimplexTree([[1,2,3], [2,3,4], [2,4,1]])
seq, n_crit = morse_seq_decroissante(K_init)

print(f"n_crit = {n_crit}")
print(f"seq = \n {seq}") 
