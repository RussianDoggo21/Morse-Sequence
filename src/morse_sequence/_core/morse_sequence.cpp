#include "morse_sequence.h"

/* ------------------------------------------------------------------------------------------------------------------------------------------------------------- */

// Constructor of the MorseSequence class
// Takes a SimplexTree as input
// O(1)
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) {
	std::cout << "MorseSequence created\n" << std::endl;
}


// Getter for the SimplexTree associated with the MorseSequence object
const SimplexTree& MorseSequence::get_simplex_tree() {
    return simplex_tree;
}

// Computes the boundary of cn with a filtration on S 
node_list MorseSequence::boundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S) {
    node_list boundary;
    simplex_t sigma = simplex_tree.full_simplex(cn); // retrieve [v₀, ..., vₚ]

    for (size_t i = 0; i < sigma.size(); ++i) {
        simplex_t face = sigma;
        face.erase(face.begin() + i); // delete the i-th vertex to create a face 

        node_ptr f = simplex_tree.find(face); // find the pointer of the face
        if (f && S.find(f) != S.end() && S.at(f)) {
            boundary.push_back(f);
        }
    }
    return boundary;
}

/*

// WAY TOO SLOW !!
node_list MorseSequence::boundary(const node_ptr& cn) {
    
    node_list boundary;
    
    faces<> bord(&simplex_tree, cn); // Iterator on the cofaces of the simplex cn

    for (auto& face : bord) { // Iterating on each face of cn
        node_ptr face_ptr = std::get<0>(face); // Retrieving the face

        if (simplex_tree.depth(face_ptr) == simplex_tree.depth(cn) + 1) { // First verification    
        boundary.push_back(face_ptr); // Adding the coface to the coboundary
        }
    }
    return boundary;
}
*/


// Computes the boundary of cn 
node_list MorseSequence::boundary(const node_ptr& cn) {
    node_list boundary;

    // if (!cn) {
    //     std::cerr << "boundary: cn is nullptr!" << std::endl;
    //     return boundary;
    // }

    simplex_t sigma = simplex_tree.full_simplex(cn); // retrieve [v₀, ..., vₚ]

    // std::cerr << "Computing boundary of simplex: ";
    // for (auto v : sigma) std::cerr << v << " ";
    // std::cerr << std::endl;

    if (sigma.size() == 1) {
        return {};
    }

    for (size_t i = 0; i < sigma.size(); ++i) {
        simplex_t face = sigma;
        face.erase(face.begin() + i); // Remove the i-th vertex

        node_ptr f = simplex_tree.find(face);
        if (f) boundary.push_back(f);
    }

    return boundary;
}




// Returns the pointers of the simplices forming the coboundary of a simplex sigma with a filtration on S
node_list MorseSequence::coboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S) {
    node_list coboundary;
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

node_list MorseSequence::coboundary(const node_ptr& cn) {
    
    node_list coboundary;

    /*
    if (!cn) {
        std::cerr << "coboundary: cn is nullptr!" << std::endl;
        return coboundary;
    }
    */
    
    cofaces<> cobord(&simplex_tree, cn); // Iterator on the cofaces of the simplex cn

    for (auto& coface : cobord) { // Iterating on each coface of cn
        node_ptr coface_ptr = std::get<0>(coface); // Retrieving the coface
        
        /*if (!coface_ptr){
            std::cerr << "coboundary : coface_ptr is nullptr!" << std::endl;
            continue;
        }
        */

        if (simplex_tree.depth(coface_ptr) == simplex_tree.depth(cn) + 1) { // First verification    
            coboundary.push_back(coface_ptr); // Adding the coface to the coboundary
        }
    }
    return coboundary;
}




// Returns the number of faces of the simplex linked to pointer cn with a filtration on S
int MorseSequence::nbboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S){
    return (this->boundary(cn,S)).size();
}

int MorseSequence::nbboundary(const node_ptr& cn){
    return (this->boundary(cn)).size();
}


// Returns the number of cofaces of the simplex linked to pointer cn with a filtration on S
int MorseSequence::nbcoboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S){
    return (this->coboundary(cn, S)).size();
}

int MorseSequence::nbcoboundary(const node_ptr& cn){
    return (this->coboundary(cn)).size();
}


// Returns the p-simplices of the simplicial complex if p is specified
// Returns all simplices if p is not specified
// Let n = simplextree.size() : O(n)
node_list MorseSequence::simplices(std::optional<int> p = std::nullopt) const {
	node_list F;

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
// Private function
node_ptr MorseSequence::find_out(const tsl::robin_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, const node_ptr& s_ptr){
    node_ptr v = nullptr;
    if (order == "decreasing"){
        for (node_ptr v0 : simplex_list){
            auto it = T.find(v0);
            if (it != T.end() && !it->second && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) + 1)) {
                v = v0;
            }       
        }
    }
    else if (order == "increasing"){
        for (node_ptr v0 : simplex_list){
            auto it = T.find(v0);
            if (it != T.end() && !it->second && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) - 1)) {
                v = v0;
            }
        }
    }
    return v;
}

// Returns the pointer of a simplex in simplex_list that satisfies a condition on T, s_ptr, and F
// Private function
node_ptr MorseSequence::find_out(const tsl::robin_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const tsl::robin_map<node_ptr, int>& F){
    node_ptr v = nullptr;
    for (node_ptr v0 : simplex_list){
        auto it = T.find(v0);
        if (it != T.end() && !it->second && F.at(v0) == F.at(s_ptr)) {
            v = v0;
        }
    }
    
    return v;
}

/*
node_list MorseSequence::get_node_list(const std::list<simplex_t>& py_list) const {
    node_list S;
    for (const simplex_t& sigma : py_list) {
        node_ptr cn = simplex_tree.find(sigma);
        if (cn == nullptr) {
            throw std::runtime_error("Simplex not found in SimplexTree");
        }
        S.push_back(cn);
    }
    return S;
}
*/


/* ------------------------------------------------------------------------------------------------------------------------------------------------------------- */

// Creation of an increasing Morse sequence
std::pair<m_sequence, int> MorseSequence::increasing(const SimplexTree& st){

    // Sort simplices by increasing dimension
    node_list K = this->simplices(); 
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return st.depth(a) < st.depth(b); 
    });

    // Declaration of variables
    m_sequence MorseSequence; // Morse sequence to be returned
    tsl::robin_map<node_ptr, bool> T; // Marks simplices used to create MorseSequence
    std::deque<node_ptr> L; // List of free tau candidates of dimension p that can form a free pair (p-1, p)
    std::unordered_map<node_ptr, int> rho; // rho[tau] == n means that tau has n faces
    int N = K.size(); // Number of simplices in st
    int i = 0; // Index to browse simplices of st
    int n_crit = 0; // Counts the number of critical simplices in MorseSequence

    T.reserve(N);
    rho.reserve(N);

    // Initialization of T and Sdict
    for (node_ptr cn : this->simplices(std::nullopt)) {
        T[cn] = false;
    }

    // Initialization of Sdict, rho, and L
    for (node_ptr cn : K) {
        int nb = this->nbboundary(cn); // number of faces of cn
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
                node_list bd = this->boundary(tau_ptr);                
                node_ptr sigma_ptr = this->find_out(T, bd, "increasing", tau_ptr);

                // Update MorseSequence and T
                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                // Update rho and L
                node_list cob1 = this->coboundary(sigma_ptr);
                node_list cob2 = this->coboundary(tau_ptr);
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
            //for (node_ptr tau_ptr : this->coboundary(sigma_ptr, Sdict)) {
            for (node_ptr tau_ptr : this->coboundary(sigma_ptr)) {
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
std::pair<m_sequence, int> MorseSequence::decreasing(const SimplexTree& st){

    // Sort simplices by decreasing dimension
    node_list K = this->simplices(); 
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return st.depth(a) > st.depth(b); 
    });

    // Declaration of variables
    m_sequence MorseSequence; // Morse sequence to be returned
    tsl::robin_map<node_ptr, bool> T; // Marks simplices used to create MorseSequence
    std::deque<node_ptr> L; // List of free sigma candidates of dimension p that can form a free pair (p, p+1)
    std::unordered_map<node_ptr, int> rho; // rho[sigma] == n means that sigma has n cofaces
    int N = K.size(); // Number of simplices in st
    int i = 0; // Index to browse simplices of st
    int n_crit = 0; // Counts the number of critical simplices in MorseSequence
 
    T.reserve(N);
    rho.reserve(N);

    // Initialization of T and Sdict 
    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
    }

    // Initialization of Sdict, rho and L
    for (node_ptr cn : K) {
        int nb = this->nbcoboundary(cn); // number of cofaces of cn
        rho[cn] = nb;
        if (nb == 1) {
            L.push_back(cn); // possible lower element of a free pair
        }
    }

    while (i < N) { // while we still haven't used all the simplices of st
        while (!L.empty()) { // while we still can add free pairs
            
            // sigma_ptr: lower element of the free pair
            node_ptr sigma_ptr = L.back(); 
            L.pop_back();
 
            
            if (rho[sigma_ptr] == 1) { // if sigma_ptr has only one coface

                // tau_ptr: upper element of the free pair
                node_list cofaces = this->coboundary(sigma_ptr);
                node_ptr tau_ptr = this->find_out(T, cofaces, "decreasing", sigma_ptr);

                // Update MorseSequence and T
                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[sigma_ptr] = true;
                T[tau_ptr] = true;

                // Update rho and L
                node_list bd1 = this->boundary(sigma_ptr);
                node_list bd2 = this->boundary(tau_ptr);
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
            for (node_ptr mu : this->boundary(sigma_ptr)) {
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
std::pair<m_sequence, int> MorseSequence::Max(const node_list& S, const tsl::robin_map<node_ptr, int>& F) {
    tsl::robin_map<node_ptr, bool> T; // Boolean dictionary: T[s] == false means s is still "available"
    tsl::robin_map<node_ptr, bool> Sdict; // Marks whether the simplex is in S
    std::deque<node_ptr> U; // Contains simplices with only one face: (sigma, tau) is a free pair
    std::unordered_map<node_ptr, int> rho; // Number of faces of a simplex
    m_sequence MorseSequence; // Result to return
    int N = S.size();
    int i = 0;
    int n_crit = 0;

    T.reserve(N);
    rho.reserve(N);
    Sdict.reserve(N);

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
                node_list boundary = this->boundary(tau_ptr, Sdict);
                node_ptr sigma_ptr = this->find_out(T, boundary, tau_ptr, F);

                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                node_list combined = this->coboundary(sigma_ptr, Sdict);
                node_list tau_coboundary = this->coboundary(tau_ptr, Sdict);
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
std::pair<m_sequence, int> MorseSequence::Min(const node_list& S, const tsl::robin_map<node_ptr, int>& F) {
    tsl::robin_map<node_ptr, bool> T; // Boolean dictionary: T[s] == false means s is still "available"
    tsl::robin_map<node_ptr, bool> Sdict;; // Marks whether the simplex is in S
    std::deque<node_ptr> U; // Contains simplices with only one coface: (sigma, tau) is a free pair
    std::unordered_map<node_ptr, int> rho; // Number of cofaces of a simplex
    m_sequence MorseSequence; // Result to return
    int N = S.size();
    int i = 0;
    int n_crit = 0;
    
    T.reserve(N);
    rho.reserve(N);
    Sdict.reserve(N);

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
                node_list cofaces = this->coboundary(sigma_ptr, Sdict);
                node_ptr tau_ptr = this->find_out(T, cofaces, sigma_ptr, F);

                MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                T[tau_ptr] = true;
                T[sigma_ptr] = true;

                node_list combined = this->boundary(sigma_ptr, Sdict);
                node_list tau_boundary = this->boundary(tau_ptr, Sdict);
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


// Print the morse_reference and the number of critical simplices
void MorseSequence::print_morse_sequence(const std::pair<m_sequence, int>& result, bool n_crit){
    
    // 1. Extract results from the function
	auto& morse_sequence = result.first;  // Morse Sequence
	int n_crit2 = result.second;  // Number of critical simplices

	// Display results (optional)
    if (n_crit == true){
        std::cout << "Number of critical points: " << n_crit2 << "\n" << std::endl;
    }
	
	for (const auto& item : morse_sequence) {
        // Check the type of the element

        if (std::holds_alternative<node_ptr>(item)) { // Critical simplex
            node_ptr face_ptr = std::get<node_ptr>(item);
            // Check here if face_ptr is not a null pointer before using it
            if (face_ptr) {
                std::cout << "Critical simplex: ";
                simplex_tree.print_simplex(std::cout, face_ptr, true);
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }
        else if (std::holds_alternative<node_pair>(item)) { // Free pair
            node_pair pair = std::get<node_pair>(item);
            
            if (pair.first && pair.second) {
                std::cout << "Pair of simplices: ";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " and ";
                simplex_tree.print_simplex(std::cout, pair.second, true);
            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }
    }  
    printf("\n"); 
}


/* 
    Computation and display of morse frames with the bitmap implementation 
*/


// Generates a map critical_simplex -> index from a Morse sequence W
node_index_map MorseSequence::generate_critical_index_map(const m_sequence& W) {
    node_index_map critical_index_map;
    size_t index = 0;
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) { // Filtration on critical simplices
            node_ptr c = std::get<node_ptr>(item);
            if (critical_index_map.find(c) == critical_index_map.end()) {
                critical_index_map[c] = index++;
            }
        }
    }
    return critical_index_map;
}

// Computes a reference map from a Morse sequence W given critical_index_map
m_frame MorseSequence::reference_map(const m_sequence& W, const node_index_map& critical_index_map) {

    m_frame ref_map;
    ref_map.reserve(W.size());

    std::size_t crit_size = critical_index_map.size();

    for (const auto& item : W) {

        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit_ptr = std::get<node_ptr>(item);

            bitmap b(crit_size);  // all zeros
            b.set(critical_index_map.at(crit_ptr));  // set the corresponding bit
            ref_map[crit_ptr] = std::move(b);
        }

        else if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);

            // ⋏(τ) = 0
            ref_map[tau_ptr] = bitmap(crit_size);

            // ∂(τ)\{σ}
            node_list bd = this->boundary(tau_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), sigma_ptr), bd.end());

            // ⋏(σ) = ⋏(∂(τ) \ {σ}) = XOR of ⋏(cn)
            bitmap acc(crit_size);  // all zeros
            for (node_ptr cn : bd) {
                acc ^= ref_map[cn];
            }
            ref_map[sigma_ptr] = std::move(acc);
        }

        else {
            printf("Error: Value of the simplices unknown");
        }
    }

    return ref_map;
}


// Computes a coreference map from a Morse sequence W given critical_index_map
m_frame MorseSequence::coreference_map(const m_sequence& W, const node_index_map& critical_index_map) {

    m_frame coref_map;
    coref_map.reserve(W.size());

    std::size_t crit_size = critical_index_map.size();

    for (auto it = W.rbegin(); it != W.rend(); ++it) {
        const auto& item = *it;

        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit_ptr = std::get<node_ptr>(item);

            bitmap b(crit_size);
            b.set(critical_index_map.at(crit_ptr));
            coref_map[crit_ptr] = std::move(b);
        }

        else if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);

            // ⋎(σ) = 0
            coref_map[sigma_ptr] = bitmap(crit_size);

            // δ(σ) \ {τ}
            node_list bd = this->coboundary(sigma_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), tau_ptr), bd.end());

            // ⋎(τ) = ⋎(δ(σ) \ {τ}) = XOR of ⋎(cn)
            bitmap acc(crit_size);
            for (node_ptr cn : bd) {
                acc ^= coref_map[cn];
            }
            coref_map[tau_ptr] = std::move(acc);
        }

        else {
            printf("Error: Value of the simplices unknown");
        }
    }

    return coref_map;
}


// Private function used in print_m_frame
void MorseSequence::print_bitmap(const bitmap& bm, const m_sequence& W, const node_index_map& critical_index_map) const {
    bool empty = true;

    for (size_t i = 0; i < bm.size(); ++i) {
        if (bm[i]) {
            empty = false;
            for (const auto& other_item : W) {
                if (std::holds_alternative<node_ptr>(other_item)) {
                    node_ptr c = std::get<node_ptr>(other_item);
                    if (critical_index_map.at(c) == i) {
                        simplex_tree.print_simplex(std::cout, c, false);
                        std::cout << " ";
                        break;
                    }
                }
            }
        }
    }

    if (empty) {
        std::cout << "{}";
    }
}



// Print a Morse Frame (reference or coreference map)
void MorseSequence::print_m_frame(const m_frame& map, const m_sequence& W, const node_index_map& critical_index_map) {
    for (const auto& item : W) {

        if (std::holds_alternative<node_ptr>(item)) { // Critical simplex
            node_ptr face_ptr = std::get<node_ptr>(item);
            if (face_ptr) {
                std::cout << "Key (Critical simplex): ";
                simplex_tree.print_simplex(std::cout, face_ptr, false);
                std::cout << " -> Value: ";
                print_bitmap(map.at(face_ptr), W, critical_index_map);
                std::cout << "\n";
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }

        else if (std::holds_alternative<node_pair>(item)) { // Free pair
            node_pair pair = std::get<node_pair>(item);

            if (pair.first && pair.second) {
                std::cout << "Pair of simplices:\n";

                std::cout << "Key (Lower pair): ";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " -> Value: ";
                print_bitmap(map.at(pair.first), W, critical_index_map);
                std::cout << "\n";

                std::cout << "Key (Upper pair): ";
                simplex_tree.print_simplex(std::cout, pair.second, false);
                std::cout << " -> Value: ";
                print_bitmap(map.at(pair.second), W, critical_index_map);
                std::cout << "\n";

            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }

        std::cout << "\n";
    }
}



/* 
    Computation and display of morse frames with the simplex implementation 
*/


// Define the symetric diffrence : A△B = AUB \ AnB 
// Private function
node_list MorseSequence::sym_diff(const node_list& A, const node_list& B){

    node_list result; // Final symmetric difference to be returned
    result.reserve(A.size() + B.size()); // To avoid re‑alloc
    std::unordered_set<node_ptr> node_buffer; // Contains the nodes seen in A 

    // Firstly, we memorize every node that is in A
    for (node_ptr cn : A){
        node_buffer.insert(cn);
    }  

    // Then we scan B
    for (node_ptr cn : B){
        if (node_buffer.erase(cn) == 0){ // If cn was not in node_buffer, it means cn is only in B              
            result.push_back(cn); // So cn is part of A△B              
        }
    }

    // Finally, whatever is still in 'node_buffer' is unique to A, the rest has been erased (case where node_buffer.erase(cn) == 1)
    result.insert(result.end(), node_buffer.begin(), node_buffer.end());

    // Remove nullptr if present in result (optional safety)
    result.erase(std::remove(result.begin(), result.end(), nullptr), result.end());

    if (result.empty()) { // If symmetric difference is empty, return list containing nullptr
        return node_list{nullptr};
    } else {
        return result; // A△B
    }
}




// Computes the reference map from a morse sequence W
m_frame0 MorseSequence::reference_map0(const m_sequence& W){
    
    m_frame0 ref_map;
    ref_map.reserve(W.size());

    // For a reference map we go through W from left to right (beginning to the end)
    for (const auto& item : W){

        // Critical simplex
        if (std::holds_alternative<node_ptr>(item)){
            node_ptr crit_ptr = std::get<node_ptr>(item); // Retrieving the simplex
            ref_map[crit_ptr] = {crit_ptr}; // ⋏(ν) = ν
        }

        // Simplex pair
        else if (std::holds_alternative<node_pair>(item)){
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item); // Retrieving the pair
            ref_map[tau_ptr] = {nullptr}; // ⋏(τ) = 0

            // Computation of the boundary ∂(τ)\{σ}
            node_list bd = this->boundary(tau_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), sigma_ptr), bd.end());

            // ⋏(σ) = ⋏(∂(τ)\{σ})
            ref_map[sigma_ptr] = {nullptr};
            for (node_ptr cn : bd){
                ref_map[sigma_ptr] = sym_diff(ref_map[sigma_ptr], ref_map[cn]);
            }
        }

        else{
            printf("Error : Value of the simplices unknown");
        } 
    }

    return ref_map;
}


// Computes the reference map from a morse sequence W
m_frame0 MorseSequence::coreference_map0(const m_sequence& W){
    
    m_frame0 coref_map;
    coref_map.reserve(W.size());

    // For a coreference map we go through W from right to left (end to start)
    for (auto it = W.rbegin(); it != W.rend(); ++it) {
        const auto& item = *it;

        // Critical simplex
        if (std::holds_alternative<node_ptr>(item)){
            node_ptr crit_ptr = std::get<node_ptr>(item); // Retrieving the simplex
            coref_map[crit_ptr] = {crit_ptr}; // ⋎(ν) = ν
        }

        // Simplex pair
        else if (std::holds_alternative<node_pair>(item)){
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item); // Retrieving the pair
            coref_map[sigma_ptr] = {nullptr}; // ⋎(σ) = 0

            // Computation of the coboundary δ(σ)\{τ}
            node_list bd = this->coboundary(sigma_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), tau_ptr), bd.end());

            //  ⋎(τ) = ⋎(δ(σ)\{τ})
            coref_map[tau_ptr] = {nullptr};
            for (node_ptr cn : bd){
                coref_map[tau_ptr] = sym_diff(coref_map[tau_ptr], coref_map[cn]);
            }
        }

        else{
            printf("Error : Value of the simplices unknown");
        } 
    }
    
    return coref_map;
}

// Print a Morse Frame (in particular a reference map or a coreference map)
void MorseSequence::print_m_frame0(m_frame0& map, const m_sequence& W){
    for (const auto& item : W){

        if (std::holds_alternative<node_ptr>(item)) { // Critical simplex
            node_ptr face_ptr = std::get<node_ptr>(item);
            // Check here if face_ptr is not a null pointer before using it
            if (face_ptr) {
                std::cout << "Key (Critical simplex): ";
                simplex_tree.print_simplex(std::cout, face_ptr, false);
                std::cout << " -> Value: ";
                for (node_ptr cn : map[face_ptr]){
                    simplex_tree.print_simplex(std::cout, cn, false);
                }  
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
            printf("\n");
        }
        else if (std::holds_alternative<node_pair>(item)) { // Free pair
            node_pair pair = std::get<node_pair>(item);
            
            if (pair.first && pair.second) {
                std::cout << "Pair of simplices: \n";
                std::cout << "Key (Lower pair)";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " -> Value: ";
                for (node_ptr cn : map[pair.first]){
                    simplex_tree.print_simplex(std::cout, cn, false);
                }  
                printf("\n");
                std::cout << "Key (Upper pair) ";
                simplex_tree.print_simplex(std::cout, pair.second, false);
                std::cout << " -> Value: ";
                for (node_ptr cn : map[pair.second]){
                    simplex_tree.print_simplex(std::cout, cn, false);
                }  
            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}