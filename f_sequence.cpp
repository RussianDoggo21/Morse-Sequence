/* IDENTIFIANT DE POINTEURS/DE SIMPLEXES : st.reindex -> st.get_vertices() ? */

#include "morse_sequence.h"
#include <deque>

// Constructeur de la classe MorseSequence
// Prend en entrée un SimplexTree
// O(1)
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) {
	std::cout << "MorseSequence créée" << std::endl;
}

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

// Renvoie le node_ptr lié au simplexe donné en paramètre
node_ptr MorseSequence::find_node(simplex_t s){
    return simplex_tree.find(s);
}

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
                node_ptr sigma_ptr;
                for (node_ptr cn : boundary){
                    if (!T[cn]){ // The condition "if not T[s]" allows us to work with a changing simplicial complex
                        sigma_ptr = cn;
                        break;
                    }
                }
                
                if (F.at(sigma_ptr) == F.at(tau_ptr)){ // Second verification on tau and sigma : is it really a free pair (by the definition of F) 
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
            node_ptr tau_ptr = U.front();
            U.pop_front();

            if (rho[tau_ptr] == 1) {
                std::vector<node_ptr> cofaces = this->coboundary(tau_ptr, Sdict);
                node_ptr sigma_ptr;
                for (node_ptr cn : cofaces) {
                    if (!T[cn]) {
                        sigma_ptr = cn;
                        break;
                    }
                }

                if (F.at(sigma_ptr) == F.at(tau_ptr)) {
                    MorseSequence.push_back(std::make_pair(tau_ptr, sigma_ptr));
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


