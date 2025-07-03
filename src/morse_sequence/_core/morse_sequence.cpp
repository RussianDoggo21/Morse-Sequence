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

// Computes the coboundary of cn with a filtration on S 
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

// Computes the coboundary of cn 
node_list MorseSequence::boundary(const node_ptr& cn) {
    node_list boundary;
    simplex_t sigma = simplex_tree.full_simplex(cn); // retrieve [v₀, ..., vₚ]
    
    if (sigma.size() == 1){
        return {};
    }

    for (size_t i = 0; i < sigma.size(); ++i) {
        simplex_t face = sigma;
        face.erase(face.begin() + i); 

        node_ptr f = simplex_tree.find(face); 
        if (f && f != nullptr) {
            boundary.push_back(f);
        }
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
    cofaces<> cobord(&simplex_tree, cn); // Iterator on the cofaces of the simplex cn

    for (auto& coface : cobord) { // Iterating on each coface of cn
        node_ptr coface_ptr = std::get<0>(coface); // Retrieving the coface
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
*/


// Computes the reference/coreference of a single simplex 
// (in particular lower simplices in reference map and upper simplices in coreference map)
// Private function
bitmap MorseSequence::Gamma(
    const node_ptr& cn, 
    const node_map& sigma2tau, 
    const node_map& tau2sigma, 
    m_frame& cache, 
    GammaMode mode, 
    const tsl::robin_map<node_ptr, std::size_t>& critical_index_map) {
    
    // Specific case : if Gamma(cn) is already cached
    auto hit = cache.find(cn);
    if (hit != cache.end()) return hit->second;

    // General case
    node_ptr sigma_ptr, tau_ptr;
    std::function<node_list(node_ptr)> traversal;
    node_ptr erase_ptr;
    node_ptr root;

    // Either we compute the reference of cn 
    if (mode == GammaMode::Reference) {

        auto it = sigma2tau.find(cn);
        auto it2 = tau2sigma.find(cn);

        // Case of critical simplex: Υ′(c) = c
        if (it == sigma2tau.end() && it2 == tau2sigma.end()) {
            bitmap bm(critical_index_map.size(), false); // We initialize bm with all bits equal to 0
            bm[critical_index_map.at(cn)] = true; // We put the bit representing cn at 1, since it's a critical simplex
            cache[cn] = bm; // We cache the result
            return bm; // We return the result
        }

        // Case of upper regular node: Υ′(τ) = 0
        // Nothing to do here
        else if (it2 != tau2sigma.end()) {
            bitmap bm(critical_index_map.size(), false); 
            cache[cn] = bm; 
            return bm;
        }
        // Case of lower regular node: Υ′(σ) = Γ(∂τ \ σ)
        else {
            sigma_ptr = cn; // The simplex on which we want to compute the reference
            tau_ptr = it->second; // The second element of the free pair (sigma_ptr, tau_ptr)
            traversal = [&](node_ptr n) { return this->boundary(n); }; // The function needed to compute the reference
            erase_ptr = sigma_ptr; // Variable generated to update cache outside of this if/else statement
            root = tau_ptr; // Variable to call traversal correctly outside of this if/else statement
        }
    } 

    // Else we compute the coreference
    else {
        auto it = tau2sigma.find(cn);
        auto it2 = sigma2tau.find(cn);

        // Case of critical simplex: Υ′′(c) = c
        if (it == tau2sigma.end() && it2 == sigma2tau.end()) {
            bitmap bm(critical_index_map.size(), false);
            bm[critical_index_map.at(cn)] = true;
            cache[cn] = bm;
            return bm;
        }
        // Case of lower regular node: Υ′′(σ) = 0
        else if (it2 != sigma2tau.end()) {
            bitmap bm(critical_index_map.size(), false);
            cache[cn] = bm;
            return bm;
        }
        // Case of upper regular node: Υ′′(τ) = Γ(∂σ \ τ)
        else {
            tau_ptr = cn; // The simplex on which we want to compute the coreference
            sigma_ptr = it->second; // The second element of the free pair (sigma_ptr, tau_ptr)
            traversal = [&](node_ptr n) { return this->coboundary(n); }; // The function needed to compute the coreference
            erase_ptr = tau_ptr; // Variable generated to update cache outside of this if/else statement
            root = sigma_ptr; // Variable to call traversal correctly outside of this if/else statement
        }
    }

    // Compute the boundary or coboundary of the root
    node_list bd = traversal(root);
    bd.erase(std::remove(bd.begin(), bd.end(), erase_ptr), bd.end());

    // Initialize result bitmap (all bits at 0)
    bitmap result(critical_index_map.size(), false);

    // Combine recursive calls with XOR
    for (node_ptr v : bd) {
        bitmap gamma_v = Gamma(v, sigma2tau, tau2sigma, cache, mode, critical_index_map);
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] ^= gamma_v[i]; // Symmetrical difference
        }
    }

    cache[erase_ptr] = result;
    return result;
}

/*
// Computes the reference/coreference of a single simplex 
// (in particular lower simplices in reference map and upper simplices in coreference map)
// Private function
node_list MorseSequence::Gamma(const node_ptr& cn, const node_map& sigma2tau, const node_map& tau2sigma, m_frame& cache, GammaMode mode) {

    // cn : simplex on which we need to compute the reference
    // sigma2tau : map containings all free pairs (sigma, tau)
    // tau2sigma : map containings all inversed free pairs (tau, sigma)
    // cache : used in the process of memoization
    // mode : to know if we compute a reference or a coreference

    //printf("cache :\n");
    //this->print_m_frame(cache);

    // Special case : gamma(cn) has already been computed before,
    // gamma(cn) is stored in cache -> We make use of memoization
    auto hit = cache.find(cn);
    if (hit != cache.end())
        return hit->second;

    // General case : we do not know gamma(cn)
    node_ptr sigma_ptr, tau_ptr;

    // Function pointer or lambda to get the traversal list (boundary or coboundary)
    std::function<node_list(node_ptr)> traversal;
    node_ptr erase_ptr;
    node_ptr root;

    // Case distinction 
    // Are we computing a reference for a lower simplex sigma_ptr or a coreference for an upper simplex 
    if (mode == GammaMode::Reference) {
        
        // Iterators over both node_map
        auto it = sigma2tau.find(cn);
        auto it2 = tau2sigma.find(cn);

        // if cn is critical
        if (it == sigma2tau.end() && it2 == tau2sigma.end()) {   
            //printf("Simplexe critique trouvé : ");
            //simplex_tree.print_simplex(std::cout, cn, true);
            cache[cn] = node_list{cn}; // σ is critical -> Γ'(σ) = {σ}
            return node_list{cn}; // We are done
        }

        // if cn is an upper simplex 
        else if (it2 != tau2sigma.end()){ 
            //printf("Simplexe supérieur trouvé : ");
            //simplex_tree.print_simplex(std::cout, cn, true);
            cache[cn] = node_list{nullptr}; // Γ'(cn) = "0"
            return node_list{nullptr}; // We are done
        }

        // if cn is a lower simplex
        else { 
            //printf("Simplexe inférieur trouvé (cas général)");
            //simplex_tree.print_simplex(std::cout, cn, true);
            sigma_ptr = cn; // The simplex on which we want to compute the reference
            tau_ptr = it->second; // The second element of the free pair (sigma_ptr, tau_ptr)
            traversal = [&](node_ptr n) { return this->boundary(n); }; // The function needed to compute the reference
            erase_ptr = sigma_ptr; // Variable generated to update cache outside of this if/else statement
            root = tau_ptr; // Variable to call traversal correctly outside of this if/else statement
        }

        

    } else { // mode == GammaMode::Coreference
        
        // iterators over both node_map
        auto it = tau2sigma.find(cn);
        auto it2 = sigma2tau.find(cn);

        // if cn is critical
        if (it == tau2sigma.end() && it2 == sigma2tau.end()) {   
            //printf("Simplexe critique trouvé : ");
            //simplex_tree.print_simplex(std::cout, cn, true);
            cache[cn] = node_list{cn}; // Γ''(σ) = {σ}
            return node_list{cn}; // We are done
        }

        // if cn is a lower simplex
        else if (it2 != sigma2tau.end()){ 
            //printf("Simplexe inférieur trouvé : ");
            //simplex_tree.print_simplex(std::cout, cn, true);
            cache[cn] = node_list{nullptr}; // Γ''(cn) = "0"
            return node_list{nullptr}; // We are done
        }

        // if cn is a upper simplex
        else{
            //printf("Simplexe supérieur trouvé (cas général)");
            //simplex_tree.print_simplex(std::cout, cn, true);
            tau_ptr = cn; // The simplex on which we want to compute the coreference
            sigma_ptr = it->second; // The second element of the free pair (sigma_ptr, tau_ptr)
            traversal = [&](node_ptr n) { return this->coboundary(n); }; // The function needed to compute the coreference
            erase_ptr = tau_ptr; // Variable generated to update cache outside of this if/else statement
            root = sigma_ptr; // Variable to call traversal correctly outside of this if/else statement

        }
    }

    // Get ∂τ\{σ} or δσ\{τ}
    node_list bd = traversal(root);
    bd.erase(std::remove(bd.begin(), bd.end(), erase_ptr), bd.end());
   
    // Initialization of the final result at "0"
    node_list result{nullptr};

    // Computation of the final result
    // Recursive call of Gamma on each node_ptr in bd to "add them up" afterwards, using the symmetrical difference
    for (node_ptr v : bd) {
        node_list gamma_v = Gamma(v, sigma2tau, tau2sigma, cache, mode); // Recursive call

        // Trivial case, gamma_v = "0" (node_list{null_ptr})
        // A △ "0" = A : No need to compute the symmetrical difference
        if (gamma_v == node_list{nullptr}) continue;

        // General case, gamma_v != "0"
        if (result == node_list{nullptr}) {
            result = gamma_v; // Special case, no need to compute the symmetrical difference
        } else {
            result = sym_diff(result, gamma_v);
        }
    }

    //printf("\n\n");

    cache[erase_ptr] = result;
    return result;
}
*/

// Generates a map critical_simplex -> index from a Morse sequence W
node_index_map MorseSequence::generate_critical_index_map(const m_sequence& W) {
    node_index_map critical_index_map;
    size_t index = 0;
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
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

    node_map sigma2tau, tau2sigma;
    m_frame ref_map;
    ref_map.reserve(W.size());

    // 1) Build sigma2tau and tau2sigma maps for free pairs
    for (const auto& item : W) {
        if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);
            sigma2tau[sigma_ptr] = tau_ptr;
            tau2sigma[tau_ptr] = sigma_ptr;
        }
    }

    // Number of critical simplices
    size_t index = critical_index_map.size();

    // 2) Initialize cache and ref_map with bitmaps for critical simplices
    m_frame cache;
    for (const auto& [c, idx] : critical_index_map) {
        bitmap b(index, false);
        b.set(idx);
        cache[c] = b;
        ref_map[c] = b;
    }

    // 3) Compute reference values for free simplices using Gamma
    for (const auto& [sigma_ptr, tau_ptr] : sigma2tau) {
        ref_map[tau_ptr] = bitmap(index, false);  // all bits set to 0
        ref_map[sigma_ptr] = Gamma(sigma_ptr, sigma2tau, tau2sigma, cache, GammaMode::Reference, critical_index_map);
    }

    return ref_map;
}

// Computes a coreference map from a Morse sequence W given critical_index_map
m_frame MorseSequence::coreference_map(const m_sequence& W, const node_index_map& critical_index_map) {

    node_map sigma2tau, tau2sigma;
    m_frame coref_map;
    coref_map.reserve(W.size());

    // 1) Build sigma2tau and tau2sigma maps for free pairs
    for (const auto& item : W) {
        if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);
            sigma2tau[sigma_ptr] = tau_ptr;
            tau2sigma[tau_ptr] = sigma_ptr;
        }
    }

    // Number of critical simplices
    size_t index = critical_index_map.size();

    // 2) Initialize cache and coref_map with bitmaps for critical simplices
    m_frame cache;
    for (const auto& [c, idx] : critical_index_map) {
        bitmap b(index, false);
        b.set(idx);
        cache[c] = b;
        coref_map[c] = b;
    }

    // 3) Compute coreference values for free simplices using Gamma
    for (const auto& [sigma_ptr, tau_ptr] : sigma2tau) {
        coref_map[sigma_ptr] = bitmap(index, false);  // all bits set to 0
        coref_map[tau_ptr] = Gamma(tau_ptr, sigma2tau, tau2sigma, cache, GammaMode::Coreference, critical_index_map);
    }

    return coref_map;
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
                // The bitmap represents which critical simplices are in the image
                const bitmap& bm = map.at(face_ptr);
                for (size_t i = 0; i < bm.size(); ++i) {
                    if (bm[i]) {
                        // Find critical simplex with index i
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
                std::cout << "\n";
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }
        else if (std::holds_alternative<node_pair>(item)) { // Free pair
            node_pair pair = std::get<node_pair>(item);

            if (pair.first && pair.second) {
                std::cout << "Pair of simplices:\n";

                // For the lower pair
                std::cout << "Key (Lower pair): ";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " -> Value: ";
                const bitmap& bm_first = map.at(pair.first);
                for (size_t i = 0; i < bm_first.size(); ++i) {
                    if (bm_first[i]) {
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
                std::cout << "\n";

                // For the upper pair
                std::cout << "Key (Upper pair): ";
                simplex_tree.print_simplex(std::cout, pair.second, false);
                std::cout << " -> Value: ";
                const bitmap& bm_second = map.at(pair.second);
                for (size_t i = 0; i < bm_second.size(); ++i) {
                    if (bm_second[i]) {
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
                std::cout << "\n";
            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }
        std::cout << "\n";
    }
}
