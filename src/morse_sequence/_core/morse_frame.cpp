#include "morse_frame.h"

// Constructor of the class
MorseFrame::MorseFrame(MorseSequence& ms) : ms(ms), simplex_tree(ms.get_simplex_tree()) {}


/* 
    Computation and display of morse frames with the bitmap implementation 
*/


// Generates a map critical_simplex -> index from a Morse sequence W
node_index_map MorseFrame::generate_critical_index_map(const m_sequence& W) {
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
m_frame MorseFrame::reference_map(const m_sequence& W, const node_index_map& critical_index_map) {

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
            node_list bd = ms.boundary(tau_ptr);
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
m_frame MorseFrame::coreference_map(const m_sequence& W, const node_index_map& critical_index_map) {

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
            node_list bd = ms.coboundary(sigma_ptr);
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
void MorseFrame::print_bitmap(const bitmap& bm, const m_sequence& W, const node_index_map& critical_index_map) const {
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
void MorseFrame::print_m_frame(const m_frame& map, const m_sequence& W, const node_index_map& critical_index_map) {
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
node_list MorseFrame::sym_diff(const node_list& A, const node_list& B){

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
m_frame0 MorseFrame::reference_map0(const m_sequence& W){
    
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
            node_list bd = ms.boundary(tau_ptr);
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
m_frame0 MorseFrame::coreference_map0(const m_sequence& W){
    
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
            node_list bd = ms.coboundary(sigma_ptr);
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
void MorseFrame::print_m_frame0(m_frame0& map, const m_sequence& W){
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