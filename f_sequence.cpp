/* IDENTIFIANT DE POINTEURS/DE SIMPLEXES : st.reindex -> st.get_vertices() ? */

#include "morse_sequence.h"

/* ------------------------------------------------------------------------------------------------------------------------------------------------------------- */

// Constructeur de la classe MorseSequence
// Prend en entrée un SimplexTree
// O(1)
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) {
	std::cout << "MorseSequence créée" << std::endl;
}

// Getter du SimplexTree lié à l'objet MorseSequence
const SimplexTree& MorseSequence::get_simplex_tree() {
    return simplex_tree;
}

// Renvoie les pointeurs des simplexes composant le bord d'un simplexe sigma avec une filtration sur S
vector<node_ptr> MorseSequence::boundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    vector<node_ptr> boundary;
    faces<> bord(&simplex_tree, cn);
    
    for (auto& face : bord) {
        node_ptr face_ptr = std::get<0>(face);

        if (simplex_tree.depth(face_ptr) == simplex_tree.depth(cn) - 1) {
            if (S.find(face_ptr) != S.end() && S.at(face_ptr)) {
                boundary.push_back(face_ptr);
            }
        }
    }
    return boundary;
}

// Renvoie les pointeurs des simplexes composant le cobord d'un simplexe sigma avec une filtration sur S
vector<node_ptr> MorseSequence::coboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S) {
    vector<node_ptr> coboundary;
    cofaces<> cobord(&simplex_tree, cn);

    for (auto& coface : cobord) {
        node_ptr coface_ptr = std::get<0>(coface);

        if (simplex_tree.depth(coface_ptr) == simplex_tree.depth(cn) + 1) {
            
            if (S.find(coface_ptr) != S.end() && S.at(coface_ptr)) {
                coboundary.push_back(coface_ptr);
            }
        }
    }
    return coboundary;
}

// Renvoie le nombre de faces du simplexe lié au pointeur cn avec une filtration sur S
int MorseSequence::nbboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S){
    return (this->boundary(cn,S)).size();
}

// Renvoie le nombre de cofaces du simplexe lié au pointeur cn avec une filtration sur S
int MorseSequence::nbcoboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S){
    return (this->coboundary(cn, S)).size();
}

// Renvoie les p-simplexes du complexe simplicial si p est précisé
// Renvoie tous les simplexes si p n'est pas précisé
// En notant n = simplextree.size() : O(n)
vector<node_ptr> MorseSequence::simplices(std::optional<int> p = std::nullopt) const {
	vector<node_ptr> F;

    // La fonction traverse requiert que lambda_f possède trois arguments de type node_ptr, idx_t et const simplex_t&
	auto lambda_f = [&F](node_ptr cn, idx_t, const simplex_t&) -> bool{
        F.push_back(cn);
        return true;
    };

	// Cas où `p` est None → Traverser tout le complexe
	if (!p) {
		auto tr = st::k_skeleton<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), simplex_tree.dimension()); // O(n) : complexité de st::k_skeleton
		traverse(tr, lambda_f); // O(n) =  O(1) * n : complexité de lambda_f * nombres de simplexes 
	}
	
	// Cas où `p` est défini → Traverser uniquement les `p`-simplexes
	else {
		auto tr = st::k_simplices<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), *p); // O(m) avec m < n
		traverse(tr, lambda_f);  // O(m) =  O(1) * m : complexité de lambda_f * nombres de simplexes 
	}

	return F;
}

// Retourne le pointeur d'un simplex qui se trouve dans simplex_list et qui vérifie une condition de T et de s_ptr
node_ptr MorseSequence::find_out(std::unordered_map<node_ptr, bool> T, std::vector<node_ptr> simplex_list, std::string order, node_ptr s_ptr){
    node_ptr v = nullptr;
    if (order == "décroissante"){
        for (node_ptr v0 : simplex_list){
            if (!T[v0] && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) + 1)) {
                v = v0;
            }
        }
    }
    else if (order == "croissante"){
        for (node_ptr v0 : simplex_list){
            if (!T[v0] && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) - 1)) {
                v = v0;
            }
        }
    }
    return v;
}

// Retourne le pointeur d'un simplex qui se trouve dans simplex_list et qui vérifie une condition de T, de s_ptr et de F
node_ptr MorseSequence::find_out(std::unordered_map<node_ptr, bool> T,std::vector<node_ptr> simplex_list, node_ptr s_ptr, const std::unordered_map<node_ptr, int>& F){
    node_ptr v = nullptr;
    for (node_ptr v0 : simplex_list){
        if (!T[v0] && F.at(v0) == F.at(s_ptr)) {
            v = v0;
        }
    }
    
    return v;
}

/* ------------------------------------------------------------------------------------------------------------------------------------------------------------- */

// Création d'une séquence de Morse croissante
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::morse_seq_crois(const SimplexTree& st)
{
    std::vector<node_ptr> K = this->simplices(); // Récupération des simplexes
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return st.depth(a) < st.depth(b); // Tri croissant par dimension
    });

    std::unordered_map<node_ptr, bool> T; // Marque les simplexes utilisés
    std::unordered_map<node_ptr, bool> Sdict; // Pour tester si un simplex est dans S
    std::deque<node_ptr> L; // Liste des candidats tau libres
    std::unordered_map<node_ptr, int> rho; // Compte les faces
    std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> MorseSequence;
    
    int N = K.size();
    int i = 0;
    int n_crit = 0;

    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        Sdict[cn] = false;
    }

    for (node_ptr cn : K) {
        Sdict[cn] = true;
        int nb = this->nbboundary(cn, Sdict);
        rho[cn] = nb;
        if (nb == 1) {
            L.push_back(cn);
        }
    }

    while (i < N) {
        while (!L.empty()) {
            node_ptr tau_ptr = L.back();
            L.pop_back();

            if (rho[tau_ptr] == 1) {
                std::vector<node_ptr> bd = this->boundary(tau_ptr, Sdict);
                node_ptr sigma_ptr = this->find_out(T, bd, "croissante", tau_ptr);

                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                std::vector<node_ptr> cob1 = this->coboundary(sigma_ptr, Sdict);
                std::vector<node_ptr> cob2 = this->coboundary(tau_ptr, Sdict);
                cob1.insert(cob1.end(), cob2.begin(), cob2.end());
                for (node_ptr mu : cob1) {
                    rho[mu] -= 1;
                    if (rho[mu] == 1) {
                        L.push_back(mu);
                    }
                }
            }
        }

        while (i < N && T[K[i]]) {
            i++;
        }

        if (i < N) {
            node_ptr sigma_ptr = K[i];
            MorseSequence.push_back(sigma_ptr);
            n_crit++;
            T[sigma_ptr] = true;
            for (node_ptr tau_ptr : this->coboundary(sigma_ptr, Sdict)) {
                rho[tau_ptr] -= 1;
                if (rho[tau_ptr] == 1) {
                    L.push_back(tau_ptr);
                }
            }
        }
    }

    return {MorseSequence, n_crit};
}


// Création d'une séquence de Morse décroissante
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::morse_seq_decrois(const SimplexTree& st){

    std::vector<node_ptr> K = this->simplices(); // Récupération des simplexes
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return st.depth(a) > st.depth(b); // Tri décroissant par dimension
    });

    std::unordered_map<node_ptr, bool> T; // Pour marquer les simplexes déjà utilisés
    std::unordered_map<node_ptr, bool> Sdict; // Pour tester l'appartenance à S
    std::deque<node_ptr> L; // Liste des paires libres potentielles
    std::unordered_map<node_ptr, int> rho; // Nombre de cofaces

    std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> MorseSequence;
    int N = K.size();
    int i = 0;
    int n_crit = 0;

    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        Sdict[cn] = false;
    }

    for (node_ptr cn : K) {
        Sdict[cn] = true;
        int nb = this->nbcoboundary(cn, Sdict); // nombre de cofaces du simplexe cn
        rho[cn] = nb;
        if (nb == 1) {
            L.push_back(cn); // sigma libre
        }
    }

    while (i < N) {
        while (!L.empty()) {
            node_ptr sigma_ptr = L.back();
            L.pop_back();
 
            
            if (rho[sigma_ptr] == 1) {
                std::vector<node_ptr> cofaces = this->coboundary(sigma_ptr, Sdict);
                node_ptr tau_ptr = this->find_out(T, cofaces, "décroissante", sigma_ptr);
                
                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[sigma_ptr] = true;
                T[tau_ptr] = true;

                std::vector<node_ptr> bd1 = this->boundary(sigma_ptr, Sdict);
                std::vector<node_ptr> bd2 = this->boundary(tau_ptr, Sdict);
                bd1.insert(bd1.end(), bd2.begin(), bd2.end());

                for (node_ptr mu : bd1) {
                    rho[mu] -= 1;
                    if (rho[mu] == 1) {
                        L.push_back(mu);
                    }
                }
                
            }
        }

        while (i < N && T[K[i]]) {
            i++;
        }

        if (i == N) break;

        node_ptr sigma_ptr = K[i];
        MorseSequence.push_back(sigma_ptr); // Simplexe critique
        n_crit++;
        T[sigma_ptr] = true;

        for (node_ptr mu : this->boundary(sigma_ptr, Sdict)) {
            rho[mu] -= 1;
            if (rho[mu] == 1) {
                L.push_back(mu);
            }
        }
    }

    return {MorseSequence, n_crit};
}



// Création d'une F-séquence maximale
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::Max(const vector<node_ptr>& S, const unordered_map<node_ptr, int>& F) {
    unordered_map<node_ptr, bool> T; // Boolean dictionnary : if T[s] == False, s is still "available to use" for the Morse Sequence
    unordered_map<node_ptr, bool> Sdict;   // Allows to check if the simplice is in S : Sdict[s] = false means it's not the case
    std::deque<node_ptr> U;  // List containing simplices v with only one coface tau : (v, tau) is a free pair for the Morse Sequence
    unordered_map<node_ptr, int> rho; // Dictionnary : if rho[s] == n, it means s has n faces
    vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> MorseSequence; // Morse Sequence to return 
    int N = S.size(); // Number of simplices in S (and st)
    int i = 0; // Incrementor to keep track of the simplices (0 <= i <= N)
    int n_crit = 0; // Counts the number of critical simplices

    for (node_ptr cn : this->simplices()){
        T[cn] = false;
        Sdict[cn] = false;
    }

    for (node_ptr cn : S){
        Sdict[cn] = true;
        int nb = this->nbboundary(cn, Sdict);
        rho[cn] = nb;
        if (nb == 1) { 
            U.push_back(cn);
        }
    }
    
    while (i<N){ //While we haven't used all simplices in S (and st)
        while (U.size() != 0){

            node_ptr tau_ptr = U.front();
            U.pop_front(); 

            if (rho[tau_ptr] == 1){                
                vector<node_ptr> boundary = this->boundary(tau_ptr, Sdict); // N.B. : boundary(tau, Sdict) should return only a single simplex 
                node_ptr sigma_ptr = this->find_out(T, boundary, tau_ptr, F);
                
                // Update of MorseSequence and T
                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                // Update of rho and then of U
                vector<node_ptr> combined = this->coboundary(sigma_ptr, Sdict);
                vector<node_ptr> tau_coboundary = this->coboundary(tau_ptr, Sdict);
                combined.insert(combined.end(), tau_coboundary.begin(), tau_coboundary.end());
                for (node_ptr mu_ptr : combined){
                    rho[mu_ptr] -= 1;
                    if (rho[mu_ptr] == 1){
                        U.push_back(mu_ptr);
                    }
                }
                
            }
        }

        // We skip all the simplices that have already been "used" in the Morse Sequence
        while (i < N && T[S[i]]){
            i += 1;
        }
        
        // Now that we have added all free pairs (loop while U), the next simplice should be a critical one
        // Update of MorseSequence, T, rho and U accordingly
        if (i<N){
            node_ptr sigma_ptr = S[i];
            MorseSequence.push_back(sigma_ptr);
            n_crit += 1;
            T[sigma_ptr] = true;
            for (node_ptr tau_ptr : this->coboundary(sigma_ptr, Sdict)){
                rho[tau_ptr] -= 1;
                if (rho[tau_ptr]==1){
                    U.push_back(tau_ptr);
                }
            }
        }
    }
    return {MorseSequence, n_crit};
}


// Création d'une F-séquence minimale
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::Min(const std::vector<node_ptr>& S, const std::unordered_map<node_ptr, int>& F) {
    std::unordered_map<node_ptr, bool> T; // Boolean dictionary: T[s] == False means s is still "available"
    std::unordered_map<node_ptr, bool> Sdict; // Marks whether the simplex is in S
    std::deque<node_ptr> U; // Contains simplices with only one face: (tau, sigma) is a free pair
    std::unordered_map<node_ptr, int> rho; // Number of cofaces of a simplex
    std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> MorseSequence; // Result to return
    int N = S.size();
    int i = 0;
    int n_crit = 0;

    // Initialization
    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        Sdict[cn] = false;
    }

    for (node_ptr cn : S) {
        Sdict[cn] = true;
        int nb = this->nbcoboundary(cn, Sdict);
        rho[cn] = nb;
        if (nb == 1) {
            U.push_back(cn);
        }
    }

    // Main loop
    while (i < N) {
        while (!U.empty()) {
            node_ptr sigma_ptr = U.front();
            U.pop_front();

            if (rho[sigma_ptr] == 1) {
                std::vector<node_ptr> cofaces = this->coboundary(sigma_ptr, Sdict);
                node_ptr tau_ptr = this->find_out(T, cofaces, sigma_ptr, F);
                

                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                std::vector<node_ptr> combined = this->boundary(sigma_ptr, Sdict);
                std::vector<node_ptr> tau_boundary = this->boundary(tau_ptr, Sdict);
                combined.insert(combined.end(), tau_boundary.begin(), tau_boundary.end());

                for (node_ptr mu_ptr : combined) {
                    rho[mu_ptr] -= 1;
                    if (rho[mu_ptr] == 1) {
                        U.push_back(mu_ptr);
                    }
                }
            }
        }

        while (i < N && T[S[i]]) {
            i++;
        }

        if (i < N) {
            node_ptr sigma_ptr = S[i];
            MorseSequence.push_back(sigma_ptr);
            n_crit++;
            T[sigma_ptr] = true;

            for (node_ptr tau_ptr : this->boundary(sigma_ptr, Sdict)) {
                rho[tau_ptr] -= 1;
                if (rho[tau_ptr] == 1) {
                    U.push_back(tau_ptr);
                }
            }
        }
    }

    // Reverse the MorseSequence before returning
    std::reverse(MorseSequence.begin(), MorseSequence.end());
    return {MorseSequence, n_crit};
}


