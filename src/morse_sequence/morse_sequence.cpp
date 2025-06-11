/* IDENTIFIANT DE POINTEURS/DE SIMPLEXES : st.reindex -> st.get_vertices() ? */

#include "morse_sequence.h"

/* ------------------------------------------------------------------------------------------------------------------------------------------------------------- */

// Constructor of the MorseSequence class
// Takes a SimplexTree as input
// O(1)
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) {
	std::cout << "MorseSequence created" << std::endl;
}

// Getter for the SimplexTree associated with the MorseSequence object
const SimplexTree& MorseSequence::get_simplex_tree() {
    return simplex_tree;
}

/* Function too costly due to the construction of the object faces<>
// Returns the pointers of the simplices forming the boundary of a simplex sigma with a filtration on S
vector<node_ptr> MorseSequence::boundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    vector<node_ptr> boundary;
    faces<> bord(&simplex_tree, cn); // Iterators on the faces of the simplex cn
    
    for (auto& face : bord) { // Iterating on each face of cn
        node_ptr face_ptr = std::get<0>(face); // Retrieving the face
        if (simplex_tree.depth(face_ptr) == simplex_tree.depth(cn) - 1) { // First verification: dim(face_ptr) == dim(cn) - 1
            if (S.find(face_ptr) != S.end() && S.at(face_ptr)) { // Second verification: face_ptr is in S and S[face_ptr] == true
                boundary.push_back(face_ptr); // Adding the face to the boundary
            }
        }
    }
    return boundary;
}
*/

vector<node_ptr> MorseSequence::boundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    vector<node_ptr> boundary;
    simplex_t sigma = simplex_tree.full_simplex(cn); // récupère [v₀, ..., vₚ]

    for (size_t i = 0; i < sigma.size(); ++i) {
        simplex_t face = sigma;
        face.erase(face.begin() + i); // enlève le i-ème sommet pour créer une face

        node_ptr f = simplex_tree.find(face); // retrouve le pointeur de la face
        if (f && S.find(f) != S.end() && S.at(f)) {
            boundary.push_back(f);
        }
    }
    return boundary;
}

/*
NE FONCTIONNE PAS
vector<node_ptr> MorseSequence::coboundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    vector<node_ptr> coboundary;

    for (const auto& child_uptr : cn->children) {
        node_ptr child = child_uptr.get(); // pointeur brut vers l'enfant

        // On vérifie que le simplexe enfant est bien dans S et actif
        if (S.find(child) != S.end() && S.at(child)) {
            coboundary.push_back(child);
        }
    }

    return coboundary;
}
*/


// Returns the pointers of the simplices forming the coboundary of a simplex sigma with a filtration on S
vector<node_ptr> MorseSequence::coboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S) {
    vector<node_ptr> coboundary;
    cofaces<> cobord(&simplex_tree, cn); // Iterator on the cofaces of the simplex cn

    for (auto& coface : cobord) { // Iterating on each coface of cn
        node_ptr coface_ptr = std::get<0>(coface); // Retrieving the coface
        if (simplex_tree.depth(coface_ptr) == simplex_tree.depth(cn) + 1) { // First verification
            if (S.find(coface_ptr) != S.end() && S.at(coface_ptr)) { // Second verification
                coboundary.push_back(coface_ptr); // Adding the coface to the coboundary
            }
        }
    }
    return coboundary;
}


// Returns the number of faces of the simplex linked to pointer cn with a filtration on S
int MorseSequence::nbboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S){
    return (this->boundary(cn,S)).size();
}

// Returns the number of cofaces of the simplex linked to pointer cn with a filtration on S
int MorseSequence::nbcoboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S){
    return (this->coboundary(cn, S)).size();
}

// Returns the p-simplices of the simplicial complex if p is specified
// Returns all simplices if p is not specified
// Let n = simplextree.size() : O(n)
vector<node_ptr> MorseSequence::simplices(std::optional<int> p = std::nullopt) const {
	vector<node_ptr> F;

    // lambda_f is required to use the function traverse() to traverse the simplicial simplex_tree
    // The parameters depth and sigma are mandatory even though they aren't used
	auto lambda_f = [&F](node_ptr cn, idx_t depth, const simplex_t& sigma) -> bool{
        F.push_back(cn);
        return true;
    };

	// Case where `p` is None → Traverse the entire complex
	if (!p) {
		auto tr = st::k_skeleton<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), simplex_tree.dimension()); // O(n): complexity of st::k_skeleton
		traverse(tr, lambda_f); // O(n) = O(1) * n: complexity of lambda_f * number of simplices 
	}
	
	// Case where `p` is defined → Traverse only the `p`-simplices
	else {
		auto tr = st::k_simplices<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), *p); // O(m) with m < n
		traverse(tr, lambda_f);  // O(m) = O(1) * m: complexity of lambda_f * number of simplices 
	}

	return F;
}

// Returns the pointer of a simplex in simplex_list that satisfies a condition on T and s_ptr
node_ptr MorseSequence::find_out(std::unordered_map<node_ptr, bool> T, std::vector<node_ptr> simplex_list, std::string order, node_ptr s_ptr){
    node_ptr v = nullptr;
    if (order == "decreasing"){
        for (node_ptr v0 : simplex_list){
            if (!T[v0] && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) + 1)) {
                v = v0;
            }
        }
    }
    else if (order == "increasing"){
        for (node_ptr v0 : simplex_list){
            if (!T[v0] && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) - 1)) {
                v = v0;
            }
        }
    }
    return v;
}

// Returns the pointer of a simplex in simplex_list that satisfies a condition on T, s_ptr, and F
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

// Creation of an increasing Morse sequence
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::increasing(const SimplexTree& st){

    // Sort simplices by increasing dimension
    std::vector<node_ptr> K = this->simplices(); 
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return st.depth(a) < st.depth(b); 
    });

    // Declaration of variables
    std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> MorseSequence; // Morse sequence to be returned
    std::unordered_map<node_ptr, bool> T; // Marks simplices used to create MorseSequence
    std::unordered_map<node_ptr, bool> Sdict; // Sdict[sigma] == True means that sigma is in S
    std::deque<node_ptr> L; // List of free tau candidates of dimension p that can form a free pair (p-1, p)
    std::unordered_map<node_ptr, int> rho; // rho[tau] == n means that tau has n faces
    int N = K.size(); // Number of simplices in st
    int i = 0; // Index to browse simplices of st
    int n_crit = 0; // Counts the number of critical simplices in MorseSequence

    // Initialization of T and Sdict
    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        Sdict[cn] = false;
    }

    // Initialization of Sdict, rho, and L
    for (node_ptr cn : K) {
        Sdict[cn] = true;
        int nb = this->nbboundary(cn, Sdict); // number of faces of cn
        rho[cn] = nb;
        if (nb == 1) {
            L.push_back(cn); // possible upper element of a free pair
        }
    }

    while (i < N) { // while we still haven't used all the simplices of st
        while (!L.empty()) { // while we can still add free pairs
            
            // tau_ptr: upper element of the free pair
            node_ptr tau_ptr = L.back(); 
            L.pop_back();

            if (rho[tau_ptr] == 1) { // if tau_ptr only has one face

                // sigma_ptr: lower element of the free pair
                std::vector<node_ptr> bd = this->boundary(tau_ptr, Sdict);
                node_ptr sigma_ptr = this->find_out(T, bd, "increasing", tau_ptr);

                // Update MorseSequence and T
                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                // Update rho and L
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

        while (i < N && T[K[i]]) { // while we still haven't used all the simplices of st
            i++;
        }

        if (i < N) { // At that point, L is empty
            node_ptr sigma_ptr = K[i]; // sigma_ptr is a critical simplex

            // Update MorseSequence, n_crit and T
            MorseSequence.push_back(sigma_ptr); 
            n_crit++;
            T[sigma_ptr] = true;

            // Update rho and L
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


// Creation of a decreasing Morse sequence
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::decreasing(const SimplexTree& st){

    // Sort simplices by decreasing dimension
    std::vector<node_ptr> K = this->simplices(); 
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return st.depth(a) > st.depth(b); 
    });

    // Declaration of variables
    std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> MorseSequence; // Morse sequence to be returned
    std::unordered_map<node_ptr, bool> T; // Marks simplices used to create MorseSequence
    std::unordered_map<node_ptr, bool> Sdict; // Sdict[sigma] == True means that sigma is in S
    std::deque<node_ptr> L; // List of free sigma candidates of dimension p that can form a free pair (p, p+1)
    std::unordered_map<node_ptr, int> rho; // rho[sigma] == n means that sigma has n cofaces
    int N = K.size(); // Number of simplices in st
    int i = 0; // Index to browse simplices of st
    int n_crit = 0; // Counts the number of critical simplices in MorseSequence
 

    // Initialization of T and Sdict 
    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        Sdict[cn] = false;
    }

    //printf("First use of coboundary...\n");

    // Initialization of Sdict, rho and L
    for (node_ptr cn : K) {
        Sdict[cn] = true;
        int nb = this->nbcoboundary(cn, Sdict); // number of cofaces of cn
        rho[cn] = nb;
        if (nb == 1) {
            L.push_back(cn); // possible lower element of a free pair
        }
    }

    //printf("Done\n");

    while (i < N) { // while we still haven't used all the simplices of st
        while (!L.empty()) { // while we still can add free pairs
            
            // sigma_ptr: lower element of the free pair
            node_ptr sigma_ptr = L.back(); 
            L.pop_back();
 
            
            if (rho[sigma_ptr] == 1) { // if sigma_ptr has only one coface
                
                // std::cout << "coboundary (free pair), i = " << i << std::endl;

                // tau_ptr: upper element of the free pair
                std::vector<node_ptr> cofaces = this->coboundary(sigma_ptr, Sdict);
                node_ptr tau_ptr = this->find_out(T, cofaces, "decreasing", sigma_ptr);

                /*
                std::cout << "sigma = ";
                st.print_simplex(std::cout, sigma_ptr, false);
                std::cout << "tau = ";
                st.print_simplex(std::cout, tau_ptr, true);
                */

                // Update MorseSequence and T
                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[sigma_ptr] = true;
                T[tau_ptr] = true;

                // Update rho and L
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

        while (i < N && T[K[i]]) { // while we still haven't used all the simplices of st
            i++;
        }

        if (i<N){ // At that point, L is empty 
            node_ptr sigma_ptr = K[i]; // sigma_ptr is a critical simplex

            // Update MorseSequence, n_crit and T
            MorseSequence.push_back(sigma_ptr); 
            n_crit++;
            T[sigma_ptr] = true;

            // Update rho and L
            for (node_ptr mu : this->boundary(sigma_ptr, Sdict)) {
                rho[mu] -= 1;
                if (rho[mu] == 1) {
                    L.push_back(mu);
                }
            }
        }
    }

    return {MorseSequence, n_crit};
}

// Build a maximal F-sequence
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::Max(const std::vector<node_ptr>& S, const std::unordered_map<node_ptr, int>& F) {
    std::unordered_map<node_ptr, bool> T; // Boolean dictionary: T[s] == false means s is still "available"
    std::unordered_map<node_ptr, bool> Sdict; // Marks whether the simplex is in S
    std::deque<node_ptr> U; // Contains simplices with only one face: (sigma, tau) is a free pair
    std::unordered_map<node_ptr, int> rho; // Number of faces of a simplex
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
        int nb = this->nbboundary(cn, Sdict);
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
                std::vector<node_ptr> boundary = this->boundary(tau_ptr, Sdict);
                node_ptr sigma_ptr = this->find_out(T, boundary, tau_ptr, F);

                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                std::vector<node_ptr> combined = this->coboundary(sigma_ptr, Sdict);
                std::vector<node_ptr> tau_coboundary = this->coboundary(tau_ptr, Sdict);
                combined.insert(combined.end(), tau_coboundary.begin(), tau_coboundary.end());

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

            for (node_ptr tau_ptr : this->coboundary(sigma_ptr, Sdict)) {
                rho[tau_ptr] -= 1;
                if (rho[tau_ptr] == 1) {
                    U.push_back(tau_ptr);
                }
            }
        }
    }

    return {MorseSequence, n_crit};
}



// Build a minimal F-sequence
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::Min(const std::vector<node_ptr>& S, const std::unordered_map<node_ptr, int>& F) {
    std::unordered_map<node_ptr, bool> T; // Boolean dictionary: T[s] == false means s is still "available"
    std::unordered_map<node_ptr, bool> Sdict; // Marks whether the simplex is in S
    std::deque<node_ptr> U; // Contains simplices with only one coface: (sigma, tau) is a free pair
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

    // Reverse the Morse sequence before returning
    std::reverse(MorseSequence.begin(), MorseSequence.end());
    return {MorseSequence, n_crit};
}
