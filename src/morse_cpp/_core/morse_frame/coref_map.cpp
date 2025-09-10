#include "morse_frame/coref_map.h"

/**
 * @brief Construct the CorefMap from the given Morse sequence.
 *
 * The process is analogous to RefMap, but:
 * - It iterates W in reverse (important for upper-to-lower consistency),
 * - It uses coboundaries instead of boundaries,
 * - The critical simplex still receives a canonical vector,
 * - For each pair (sigma, tau), sigma gets the zero vector and tau receives an XOR over its coboundary,
 *   excluding the one paired with sigma.
 *
 * @param ms The Morse sequence context.
 * @param W The Morse sequence.
 */
CorefMap::CorefMap(MorseSequence& ms, m_sequence W) : Utils(ms, W), W(W) {
    for (auto it = this->W.rbegin(); it != this->W.rend(); ++it) {
        const auto& item = *it;
        FullBitArrayElement element;
        if (std::holds_alternative<node_ptr>(item)) {
            // Critical simplex: assign canonical unit vector
            node_ptr crit_ptr = std::get<node_ptr>(item);
            bitmap b(dim_crit);
            b.set(critToIndex.at(crit_ptr));
            add(crit_ptr, b);
            element.is_critical = true;
            element.critical_simplex[crit_ptr] = b;
        }

        else if (std::holds_alternative<node_pair>(item)) {
            // Paired simplices (sigma, tau) — handle in reverse
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);

            // Sigma (lower) gets the zero vector
            add(sigma_ptr, bitmap(dim_crit));

            // Compute XOR of all coboundary cofaces of sigma, excluding tau
            node_list cbd = ms.coboundary(sigma_ptr);
            cbd.erase(std::remove(cbd.begin(), cbd.end(), tau_ptr), cbd.end());

            bitmap b_tau(dim_crit);
            for (node_ptr cn : cbd) {
                b_tau ^= get(cn);  // XOR with known values
            }
            // Tau (upper) gets the result
            add(tau_ptr, b_tau);
            element.is_critical = false;
            element.pair = {{sigma_ptr, bitmap(dim_crit)}, {tau_ptr, b_tau}};
        }
        else {
            throw std::runtime_error("Unknown simplex type encountered in Morse sequence.");
        }
        this->full_bitarray.push_back(element);
    }
}

/**
 * @brief Computes the copersistence pairs and essential simplices from the coreference map.
 * 
 * The copersistence is a list of pairs ((ν, σ), λ) where ν is a critical simplex,
 * σ is the critical simplex paired with ν, and λ is the persistence value (difference in weights).
 * Essential simplices are those with infinite persistence.
 * 
 * @return A pair consisting of:
 *         - A list of essential critical simplices (node_list).
 *         - A vector of persistence pairs (vector of pairs ((node_ptr, node_ptr), int)) sorted by persistence ascending.
 * 
 * @note The copersistence result is analogous to the persistence result computed from the reference map
 *       when both maps are constructed from the same Morse sequence.
 */

// To test on a copy of the original corefmap
// CorefMap morse_frame = coref_map;
// morse_frame.copersistence();
std::pair<node_list, std::vector<std::pair<node_pair, int>>> CorefMap::copersistence() {
    node_list essential;
    std::vector<std::pair<node_pair, int>> pairs;

    // Iterate over critical simplices in reverse order
    for (auto it = critics.rbegin(); it != critics.rend(); ++it) {
        node_ptr sigma = *it;

        // Compute the coboundary vector of sigma
        node_list coboundary = ms.coboundary(sigma);
        bitmap bsigma(dim_crit);
        for (node_ptr tau : coboundary) {
            bsigma ^= get(tau);  // XOR the bitarrays from the coboundary
        }

        if (bsigma.any()) {
            // Find the critical simplex nu that minimizes F[tau]
            double min_value = std::numeric_limits<double>::infinity();
            node_ptr nu;
            node_stack F = ms.get_stack();
            for (size_t i = bsigma.find_first(); i != boost::dynamic_bitset<>::npos; i = bsigma.find_next(i)) {
                node_ptr tau = indexToCrit[i];
                if (F[tau] <= min_value) {
                    min_value = F[tau];
                    nu = tau;
                }
            }

            // Store the copersistence pair ((sigma, nu), F[nu] - F[sigma])
            pairs.push_back({{sigma, nu}, F[nu] - F[sigma]});

            // Update all bitarrays to account for the copersistence pair
            int index_sigma = critToIndex[sigma];
            int index_nu = critToIndex[nu];
            transform_bitarrays([&](const bitmap& ba) -> bitmap {
                return update_bitarray(ba, index_sigma, index_nu, bsigma);
            });
        }
    }

    // Sort the copersistence pairs by increasing persistence
    std::sort(pairs.begin(), pairs.end(),
        [](const std::pair<node_pair, int>& a, const std::pair<node_pair, int>& b) {
            return a.second < b.second;
        });

    // Identify essential critical simplices
    for (node_ptr sigma : critics) {
        int idx = critToIndex[sigma];
        if (get(sigma)[idx]) {  // If sigma is still present in its bitarray
            essential.push_back(sigma);
        }
    }

    return {essential, pairs};
}

/*
tsl::robin_map<node_ptr, bitmap> CorefMap::get_full_bitarray() {
    tsl::robin_map<node_ptr, bitmap> full_bitarray;
    tsl::robin_map<node_ptr, bitmap> bitarray = get_bitarray();

    //std::cout << "DEBUG: Building full_bitarray..." << std::endl;

    std::cout << "DEBUG: W size = " << this->W.size() << std::endl;
    for (const auto& item : this->W) {
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit_ptr = std::get<node_ptr>(item);
            node_ptr root = _find(crit_ptr);
          
            full_bitarray[crit_ptr] = bitarray.at(root);
        }
        else if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);
            node_ptr sigma_root = _find(sigma_ptr);
            node_ptr tau_root = _find(tau_ptr);
           
            full_bitarray[sigma_ptr] = bitarray.at(sigma_root);
            full_bitarray[tau_ptr] = bitarray.at(tau_root);
        }
    }

    return full_bitarray;
}


void CorefMap::print_full_bitarray(const m_sequence& W) {
    // Récupère la m_frame complète (racines + paires)
    tsl::robin_map<node_ptr, bitmap> full_bitarray = get_full_bitarray();

    // Affiche chaque entrée de la full_bitarray
    for (const auto& [key_ptr, bm] : full_bitarray) {
        if (!key_ptr) {
            std::cout << "Null key encountered!" << std::endl;
            continue;
        }

        // Affiche la clé (simplex)
        std::cout << "Key: ";
        simplex_tree.print_simplex(std::cout, key_ptr, false);
        std::cout << " -> Value: ";

        // Affiche le bitmap associé
        print_bitmap(bm, W);
        std::cout << "\n";
    }
}

*/
