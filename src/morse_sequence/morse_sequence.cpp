#include "morse_sequence.h"
using m_sequence = std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>;
using node_list = std::vector<node_ptr>;
using morse_frame = std::unordered_map<node_ptr, node_list>;
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

/*
//Function too costly due to the construction of the object faces<>
// Returns the pointers of the simplices forming the boundary of a simplex sigma with a filtration on S
node_list MorseSequence::boundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    node_list boundary;
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


node_list MorseSequence::boundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    node_list boundary;
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
    

node_list MorseSequence::boundary(node_ptr cn) {
    node_list boundary;
    simplex_t sigma = simplex_tree.full_simplex(cn); // récupère [v₀, ..., vₚ]
    
    if (sigma.size() == 1){
        return {};
    }

    for (size_t i = 0; i < sigma.size(); ++i) {
        simplex_t face = sigma;
        face.erase(face.begin() + i); // enlève le i-ème sommet pour créer une face

        node_ptr f = simplex_tree.find(face); // retrouve le pointeur de la face
        if (f && f != nullptr) {
            boundary.push_back(f);
        }
    }
    return boundary;
}

/*
NE FONCTIONNE PAS
node_list MorseSequence::coboundary(node_ptr cn, const std::unordered_map<node_ptr, bool>& S) {
    node_list coboundary;

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
node_list MorseSequence::coboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S) {
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

node_list MorseSequence::coboundary(node_ptr cn) {
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
int MorseSequence::nbboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S){
    return (this->boundary(cn,S)).size();
}

int MorseSequence::nbboundary(node_ptr cn){
    return (this->boundary(cn)).size();
}

// Returns the number of cofaces of the simplex linked to pointer cn with a filtration on S
int MorseSequence::nbcoboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S){
    return (this->coboundary(cn, S)).size();
}

int MorseSequence::nbcoboundary(node_ptr cn){
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
node_ptr MorseSequence::find_out(const std::unordered_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, node_ptr s_ptr){
    node_ptr v = nullptr;
    if (order == "decreasing"){
        for (node_ptr v0 : simplex_list){
            auto it = T.find(v0);
            //if (!T[v0] && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) + 1)) {
            if (it != T.end() && !it->second && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) + 1)) {
                v = v0;
            }       
        }
    }
    else if (order == "increasing"){
        for (node_ptr v0 : simplex_list){
            auto it = T.find(v0);
            //if (!T[v0] && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) - 1)) {
             if (it != T.end() && !it->second && simplex_tree.depth(v0) == (simplex_tree.depth(s_ptr) - 1)) {
                v = v0;
            }
        }
    }
    return v;
}

// Returns the pointer of a simplex in simplex_list that satisfies a condition on T, s_ptr, and F
node_ptr MorseSequence::find_out(const std::unordered_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const std::unordered_map<node_ptr, int>& F){
    node_ptr v = nullptr;
    for (node_ptr v0 : simplex_list){
        auto it = T.find(v0);
        //if (!T[v0] && F.at(v0) == F.at(s_ptr)) {
        if (it != T.end() && !it->second && F.at(v0) == F.at(s_ptr)) {
            v = v0;
        }
    }
    
    return v;
}

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
    std::unordered_map<node_ptr, bool> T; // Marks simplices used to create MorseSequence
    //std::unordered_map<node_ptr, bool> Sdict; // Sdict[sigma] == True means that sigma is in S
    std::deque<node_ptr> L; // List of free tau candidates of dimension p that can form a free pair (p-1, p)
    std::unordered_map<node_ptr, int> rho; // rho[tau] == n means that tau has n faces
    int N = K.size(); // Number of simplices in st
    int i = 0; // Index to browse simplices of st
    int n_crit = 0; // Counts the number of critical simplices in MorseSequence

    T.reserve(N);
    rho.reserve(N);
    //Sdict.reserve(N);

    // Initialization of T and Sdict
    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        //Sdict[cn] = false;
    }

    // Initialization of Sdict, rho, and L
    for (node_ptr cn : K) {
        //Sdict[cn] = true;
        //int nb = this->nbboundary(cn, Sdict); // number of faces of cn
        int nb = this->nbboundary(cn); // number of faces of cn
        rho[cn] = nb;
        if (nb == 1) {
            L.push_back(cn); // possible upper element of a free pair
        }
    }

    /*
    for (node_ptr cn : K){
        printf("nbboundary of ");
        st.print_simplex(std::cout, cn, false);
        std::cout << " = " << rho[cn];
        printf(" Boundary of ");
        st.print_simplex(std::cout, cn, false);
        printf(": ");
        for (node_ptr c : this->boundary(cn)){
            st.print_simplex(std::cout, c, false);
            if (c == nullptr){
                printf(" null ptr detected in boundary of ");
                st.print_simplex(std::cout, cn, false);
            }
        }
        printf("\n");
    }
    */

    //printf("increasing : Initialisation terminée\n ");
    while (i < N) { // while we still haven't used all the simplices of st
        //std::cout << "i = " << i << std::endl;
        //std::cout << "size of L = " << L.size() << std::endl;
        try {
            while (!L.empty()) { // while we can still add free pairs
                //printf("Entrée dans L");
                // tau_ptr: upper element of the free pair
                node_ptr tau_ptr = L.back(); 
                L.pop_back();

                /*
                printf("\ntau = ");
                st.print_simplex(std::cout, tau_ptr, true);
                printf("rho[");
                st.print_simplex(std::cout, tau_ptr, false);
                std::cout << "] = " << rho[tau_ptr] <<  " \n";
                */

                if (rho[tau_ptr] == 1) { // if tau_ptr only has one face

                    // sigma_ptr: lower element of the free pair
                    
                    //printf("Boundary");
                    //node_list bd = this->boundary(tau_ptr, Sdict);
                    node_list bd = this->boundary(tau_ptr);
                    
                    /*
                    printf("Boundary of ");
                    st.print_simplex(std::cout, tau_ptr, false);
                    printf(": ");
                    for (node_ptr cn : bd){
                        st.print_simplex(std::cout, cn, false);
                        std::cout << "T value = " << T[cn] << "\n";
                    }
                    */
                    

                    node_ptr sigma_ptr = this->find_out(T, bd, "increasing", tau_ptr);

                    
                    //printf("sigma = ");
                    //st.print_simplex(std::cout, sigma_ptr, true);
                    

                    // Update MorseSequence and T
                    MorseSequence.push_back(std::make_pair(sigma_ptr, tau_ptr));
                    T[tau_ptr] = true;
                    T[sigma_ptr] = true;

                    // Update rho and L
                    //printf("Update rho and L");
                    //std::node_list cob1 = this->coboundary(sigma_ptr, Sdict);
                    //std::node_list cob2 = this->coboundary(tau_ptr, Sdict);
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
                //printf("Update rho and L");
                //for (node_ptr tau_ptr : this->coboundary(sigma_ptr, Sdict)) {
                for (node_ptr tau_ptr : this->coboundary(sigma_ptr)) {
                    rho[tau_ptr] -= 1;
                    if (rho[tau_ptr] == 1) {
                        L.push_back(tau_ptr);
                    }
                }
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Exception attrapée dans L : " << e.what() << std::endl;
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
    std::unordered_map<node_ptr, bool> T; // Marks simplices used to create MorseSequence
    //std::unordered_map<node_ptr, bool> Sdict; // Sdict[sigma] == True means that sigma is in S
    std::deque<node_ptr> L; // List of free sigma candidates of dimension p that can form a free pair (p, p+1)
    std::unordered_map<node_ptr, int> rho; // rho[sigma] == n means that sigma has n cofaces
    int N = K.size(); // Number of simplices in st
    int i = 0; // Index to browse simplices of st
    int n_crit = 0; // Counts the number of critical simplices in MorseSequence
 
    T.reserve(N);
    rho.reserve(N);
    //Sdict.reserve(N);

    // Initialization of T and Sdict 
    for (node_ptr cn : this->simplices()) {
        T[cn] = false;
        //Sdict[cn] = false;
    }

    //printf("First use of coboundary...\n");

    // Initialization of Sdict, rho and L
    for (node_ptr cn : K) {
        //Sdict[cn] = true;
        //int nb = this->nbcoboundary(cn, Sdict); // number of cofaces of cn
        int nb = this->nbcoboundary(cn); // number of cofaces of cn
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
                //std::node_list cofaces = this->coboundary(sigma_ptr, Sdict);
                node_list cofaces = this->coboundary(sigma_ptr);
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
                //node_list bd1 = this->boundary(sigma_ptr, Sdict);
                //node_list bd2 = this->boundary(tau_ptr, Sdict);
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
            //for (node_ptr mu : this->boundary(sigma_ptr, Sdict)) {
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
std::pair<m_sequence, int> MorseSequence::Max(const node_list& S, const std::unordered_map<node_ptr, int>& F) {
    std::unordered_map<node_ptr, bool> T; // Boolean dictionary: T[s] == false means s is still "available"
    std::unordered_map<node_ptr, bool> Sdict; // Marks whether the simplex is in S
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
std::pair<m_sequence, int> MorseSequence::Min(const node_list& S, const std::unordered_map<node_ptr, int>& F) {
    std::unordered_map<node_ptr, bool> T; // Boolean dictionary: T[s] == false means s is still "available"
    std::unordered_map<node_ptr, bool> Sdict; // Marks whether the simplex is in S
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
void MorseSequence::print_morse_sequence(std::pair<m_sequence, int> result, bool n_crit){
    
    // 1. Extract results from the function
	auto& morse_sequence = result.first;  // Morse Sequence
	int n_crit2 = result.second;  // Number of critical simplices

	// Display results (optional)
    if (n_crit == true){
        std::cout << "Number of critical points: " << n_crit2 << "\n" << std::endl;
    }
	
	for (const auto& item : morse_sequence) {
        // Check the type of the element
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr face_ptr = std::get<node_ptr>(item);
            // Check here if face_ptr is not a null pointer before using it
            if (face_ptr) {
                std::cout << "Critical simplex: ";
                simplex_tree.print_simplex(std::cout, face_ptr, true);
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }
        else if (std::holds_alternative<std::pair<node_ptr, node_ptr>>(item)) {
            std::pair<node_ptr, node_ptr> pair = std::get<std::pair<node_ptr, node_ptr>>(item);
            
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


// Define the symetric diffrence : A△B = AUB \ AnB
node_list sym_diff(const node_list& A, const node_list& B){

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

    // If symmetric difference is empty, return list containing nullptr
    if (result.empty()) {
        return node_list{nullptr};
    } else {
        return result; // A△B
    }
}

// Computes a reference map from a morse sequence
morse_frame MorseSequence::reference_map(const m_sequence& W){

    using pair_type = std::pair<node_ptr, node_ptr>;
    morse_frame ref_map;
    ref_map.reserve(W.size());           

    /*----- 1. critical simplices first ------------------------------------*/
    for (const auto& item : W){
        if (std::holds_alternative<node_ptr>(item)){
            node_ptr cn = std::get<node_ptr>(item);
            ref_map[cn] = {cn};            // f(cn) = cn
        }
    }

    /*----- 2. free pairs after ----------------------------------------------*/
    for (const auto& item : W){
        if (std::holds_alternative<pair_type>(item)){

            auto [sigma_ptr, tau_ptr] = std::get<pair_type>(item);

            ref_map[tau_ptr] = node_list{nullptr};  // ref_map[tau_ptr] = 0

            // Compute ∂τ\{σ}
            node_list bd = boundary(tau_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), sigma_ptr), bd.end());
            
            // Initialize ref_map[sigma_ptr] = 0
            ref_map[sigma_ptr] = node_list{nullptr};  

            for (node_ptr cn : bd) {
                auto it = ref_map.find(cn);
                if (it == ref_map.end() || it->second == node_list{nullptr}) {
                    continue;  // skip undefined or sentinel cases
                }
                ref_map[sigma_ptr] = sym_diff(ref_map[sigma_ptr], it->second);
            }
        }
    }
    return ref_map;
}


// Computes a coreference map from a morse sequence
morse_frame MorseSequence::coreference_map(const m_sequence& W){

    using pair_type = std::pair<node_ptr, node_ptr>;
    morse_frame ref_map;
    ref_map.reserve(W.size());           

    /*----- 1. critical simplices first ------------------------------------*/
    for (const auto& item : W){
        if (std::holds_alternative<node_ptr>(item)){
            node_ptr cn = std::get<node_ptr>(item);
            ref_map[cn] = {cn};            // f(cn) = cn
        }
    }

    /*----- 2. free pairs after ----------------------------------------------*/
    for (const auto& item : W){
        if (std::holds_alternative<pair_type>(item)){

            auto [sigma_ptr, tau_ptr] = std::get<pair_type>(item);

            ref_map[sigma_ptr] = node_list{nullptr};  // ref_map[sigma_ptr] = 0

            // Compute ∂τ\{σ}
            node_list cbd = coboundary(sigma_ptr);
            cbd.erase(std::remove(cbd.begin(), cbd.end(), tau_ptr), cbd.end());
            
            // Initialize ref_map[tau_ptr] = 0
            ref_map[tau_ptr] = node_list{nullptr};  

            for (node_ptr cn : cbd) {
                auto it = ref_map.find(cn);
                if (it == ref_map.end() || it->second == node_list{nullptr}) {
                    continue;  // skip undefined or sentinel cases
                }
                ref_map[tau_ptr] = sym_diff(ref_map[tau_ptr], it->second);
            }
        }
    }
    return ref_map;
}

// Print a Morse Frame (in particular a reference map or a coreference map)
void MorseSequence::print_morse_frame(morse_frame map){
    for (const auto& [key, val] : map) {
        std::cout << "  Key: ";
        simplex_tree.print_simplex(std::cout, key, false);
        std::cout << " -> Value: ";
        for (node_ptr cn : val){
            simplex_tree.print_simplex(std::cout, cn, false);
        }  
        std::cout << "\n";
    }
    printf("\n");
}


