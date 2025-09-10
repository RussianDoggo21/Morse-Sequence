#include "morse_sequence/morse_sequence.h"

/**
 * @brief Construct a MorseSequence from a SimplexTree.
 * @param st The simplex tree to use.
 */
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) {
	std::cout << "MorseSequence created\n" << std::endl;
    for (node_ptr cn: this->simplices(std::nullopt) ){
        F[cn] = 0;
    }
}


/**
 * @brief Get the underlying simplex tree.
 * @return A const reference to the simplex tree.
 */
const SimplexTree& MorseSequence::get_simplex_tree() {
    return simplex_tree;
}

/**
 * @brief Compute the boundary of a node within a subset.
 * @param cn The node to compute the boundary for.
 * @param S A map indicating which nodes are in the subset.
 * @return A list of boundary nodes in the subset.
 */
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


/**
 * @brief Compute the full boundary of a node.
 * @param cn The node to compute the boundary for.
 * @return A list of boundary nodes.
 */
node_list MorseSequence::boundary(const node_ptr& cn) {
    node_list boundary;

    simplex_t sigma = simplex_tree.full_simplex(cn); // retrieve [v₀, ..., vₚ]

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

/**
 * @brief Compute the coboundary of a node within a subset.
 * @param cn The node to compute the coboundary for.
 * @param S A map indicating which nodes are in the subset.
 * @return A list of coboundary nodes in the subset.
 */
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

/**
 * @brief Compute the full coboundary of a node.
 * @param cn The node to compute the coboundary for.
 * @return A list of coboundary nodes.
 */
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

/**
 * @brief Count the number of boundary nodes in the subset.
 * @param cn The node whose boundary is considered.
 * @param S A map of nodes in the subset.
 * @return The number of boundary nodes in the subset.
 */
int MorseSequence::nbboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S){
    return (this->boundary(cn,S)).size();
}

/**
 * @brief Get the full number of boundary nodes.
 * @param cn The node to inspect.
 * @return The size of the boundary list.
 */
int MorseSequence::nbboundary(const node_ptr& cn){
    return (this->boundary(cn)).size();
}


/**
 * @brief Count the number of coboundary nodes in the subset.
 * @param cn The node whose coboundary is considered.
 * @param S A map of nodes in the subset.
 * @return The number of coboundary nodes in the subset.
 */
int MorseSequence::nbcoboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S){
    return (this->coboundary(cn, S)).size();
}

/**
 * @brief Get the full number of coboundary nodes.
 * @param cn The node to inspect.
 * @return The size of the coboundary list.
 */
int MorseSequence::nbcoboundary(const node_ptr& cn){
    return (this->coboundary(cn)).size();
}


/**
 * @brief Get all simplices of the optional given dimension.
 * @param p The optional dimension to filter simplices.
 * @return A list of simplices.
 */
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


/**
 * @brief Find a simplex pointer in simplex_list that satisfies a condition based on T and s_ptr.
 * 
 * This function searches through simplex_list for a simplex v0 such that:
 * - v0 is in map T with a false boolean value (i.e., !T[v0]),
 * - and the depth condition holds:
 *   * if order == "decreasing": depth(v0) == depth(s_ptr) + 1
 *   * if order == "increasing": depth(v0) == depth(s_ptr) - 1
 * 
 * The function returns the last matching simplex pointer found.
 * 
 * @param T A map from simplex pointers to boolean flags.
 * @param simplex_list The list of simplexes to search.
 * @param order A string specifying the depth condition ("increasing" or "decreasing").
 * @param s_ptr The reference simplex pointer used to compare depths.
 * @return node_ptr The pointer to the found simplex, or nullptr if none matches.
 */
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

/**
 * @brief Find a simplex pointer in simplex_list that satisfies a condition based on T, s_ptr, and F.
 * 
 * This function searches through simplex_list for a simplex v0 such that:
 * - v0 is in map T with a false boolean value,
 * - and F[v0] == F[s_ptr].
 * 
 * The function returns the last matching simplex pointer found.
 * 
 * @param T A map from simplex pointers to boolean flags.
 * @param simplex_list The list of simplexes to search.
 * @param s_ptr The reference simplex pointer used to compare values in F.
 * @param F A map from simplex pointers to integer values.
 * @return node_ptr The pointer to the found simplex, or nullptr if none matches.
 */
node_ptr MorseSequence::find_out(const tsl::robin_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const node_stack& F){
    node_ptr v = nullptr;
    for (node_ptr v0 : simplex_list){
        auto it = T.find(v0);
        if (it != T.end() && !it->second && F.at(v0) == F.at(s_ptr)) {
            v = v0;
        }
    }
    
    return v;
}

/**
 * @brief Get the node_stack values associated to the simplicial complex.
 *
 * Returns a constant reference to the internal map that associates
 * each simplex (represented by a node pointer) to its filtration value
 * or level in the Morse sequence.
 *
 * This data is typically used during the computation of persistence pairs
 * or essential simplices.
 *
 * @return A constant reference to the map from simplex nodes to integer values.
 */
const node_stack& MorseSequence::get_stack() const{
        return F;
}

/**
 * @brief Update the internal node_stack map with a new mapping.
 * 
 * This function replaces the current node_stack (map associating nodes to integers)
 * with the provided new node_stack.
 * 
 * @param new_F The new node_stack mapping to set. 
 */
void MorseSequence::update_stack(node_stack new_F){
    F = new_F;
}


/* ------------------------------------------------------------------------------------------------------------------------------------------------------------- */

/**
 * @brief Constructs an increasing discrete Morse sequence from a simplex tree.
 * 
 * This method computes a discrete Morse sequence where the matching is done by
 * pairing simplices with one face, going from lower to higher dimension (increasing).
 * It uses a queue of candidate simplices that can form free pairs and iteratively
 * constructs the Morse sequence until all simplices are processed.
 * 
 * Algorithm overview:
 * - Sort simplices by increasing dimension.
 * - Initialize data structures for bookkeeping: T (used simplices), rho (number of faces), L (free pairs candidates).
 * - While there remain unused simplices, try to add free pairs using simplices with exactly one face.
 * - If no free pairs can be added, mark the next unused simplex as critical.
 * 
 * @param st The SimplexTree from which the Morse sequence is derived.
 * @return std::pair<m_sequence, int> A pair consisting of:
 *   - the Morse sequence as a vector of pairs (free pairs) or single simplices (critical simplices),
 *   - the count of critical simplices.
 */
std::pair<m_sequence, int> MorseSequence::increasing(){

    // Sort simplices by increasing dimension
    node_list K = this->simplices(); 
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return simplex_tree.depth(a) < simplex_tree.depth(b); 
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


/**
 * @brief Constructs a decreasing discrete Morse sequence from a simplex tree.
 * 
 * This method computes a discrete Morse sequence where the matching is done by
 * pairing simplices with one coface, going from higher to lower dimension (decreasing).
 * The algorithm is symmetric to the increasing case but uses cofaces instead of faces,
 * and simplices are sorted by decreasing dimension.
 * 
 * Algorithm overview:
 * - Sort simplices by decreasing dimension.
 * - Initialize bookkeeping structures: T, rho (number of cofaces), L.
 * - Iteratively add free pairs where simplices have exactly one coface.
 * - Mark simplices as critical if no free pairs are possible.
 * 
 * @param st The SimplexTree from which the Morse sequence is derived.
 * @return std::pair<m_sequence, int> A pair consisting of:
 *   - the Morse sequence as a vector of pairs (free pairs) or single simplices (critical simplices),
 *   - the count of critical simplices.
 */
std::pair<m_sequence, int> MorseSequence::decreasing(){

    // Sort simplices by decreasing dimension
    node_list K = this->simplices(); 
    std::sort(K.begin(), K.end(), [&](const node_ptr& a, const node_ptr& b){
        return simplex_tree.depth(a) > simplex_tree.depth(b); 
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


/**
 * @brief Build a maximal F-sequence on a given subset of simplices.
 *
 * This function constructs a maximal Morse sequence on the simplices in the cosimplicial complex S,
 * using the values in the dictionary F to guide the pairing of simplices.
 *
 * @param S A cosimplicial complex.
 * @return A pair consisting of:
 *   - The maximal Morse sequence (vector of pairs of simplices or critical simplices).
 *   - The number of critical simplices found.
 */
//std::pair<m_sequence, int> MorseSequence::Max(const node_list& S, const node_stack& F) {
std::pair<m_sequence, int> MorseSequence::Max(const node_list& S) {
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


/**
 * @brief Build a minimal F-sequence on a given subset of simplices.
 *
 * This function constructs a minimal Morse sequence on the simplices in the cosimplicial complex S,
 * using the values in the dictionary F to guide the pairing of simplices.
 *
 * @param S A cosimplicial complex.
 * @return A pair consisting of:
 *   - The minimal Morse sequence (vector of pairs of simplices or critical simplices).
 *   - The number of critical simplices found.
 */
//std::pair<m_sequence, int> MorseSequence::Min(const node_list& S, const node_stack& F) {
std::pair<m_sequence, int> MorseSequence::Min(const node_list& S) {
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

/**
 * @brief Print the Morse sequence in a readable format.
 *
 * This function iterates over each element in the Morse sequence and prints it.
 * It distinguishes between critical simplices and free pairs, printing each in a specific format.
 *
 * @param morse_sequence The Morse sequence to print, which can contain both critical simplices
 *                        (as node_ptr) and free pairs (as node_pair).
 */
void MorseSequence::print_morse_sequence0(m_sequence morse_sequence){

    // Iterate over each element in the Morse sequence
    for (const auto& item : morse_sequence) {
        // If the item is a critical simplex (single node_ptr)
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr simplex = std::get<node_ptr>(item);
            if (simplex) {
                std::cout << "Critical simplex: ";
                simplex_tree.print_simplex(std::cout, simplex, true);
            } else {
                std::cout << "Null pointer encountered in critical simplex!" << std::endl;
            }
        }
        // If the item is a free pair (pair of node_ptr)
        else if (std::holds_alternative<node_pair>(item)) {
            node_pair pair = std::get<node_pair>(item);
            if (pair.first && pair.second) {
                std::cout << "Pair of simplices: ";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " and ";
                simplex_tree.print_simplex(std::cout, pair.second, true);
            } else {
                std::cout << "Null pointer encountered in simplex pair!" << std::endl;
            }
        }
    }

    std::cout << std::endl;
}


/**
 * @brief Print the Morse sequence and optionally the number of critical simplices.
 *
 * This function takes the result of a Morse sequence computation and prints its content
 * in a readable format. It uses `print_morse_sequence0` to handle the actual printing of the sequence.
 * Critical simplices and free pairs are handled separately.
 *
 * @param result The Morse sequence result: a pair consisting of
 *               - first: the sequence of simplices or pairs,
 *               - second: the number of critical simplices.
 * @param print_crit If true, print the number of critical simplices.
 */
void MorseSequence::print_morse_sequence(const std::pair<m_sequence, int>& result, bool print_crit) {
    // Extract the Morse sequence and critical simplex count
    const auto& morse_sequence = result.first;
    int n_crit = result.second;

    // Optionally print the number of critical simplices
    if (print_crit) {
        std::cout << "Number of critical points: " << n_crit << "\n" << std::endl;
    }

    this->print_morse_sequence0(morse_sequence);
}
