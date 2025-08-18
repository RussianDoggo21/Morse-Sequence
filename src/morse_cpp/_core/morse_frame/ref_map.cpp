#include "morse_frame/ref_map.h"

/**
 * @brief Construct a RefMap object and compute the Morse frame.
 *
 * Each critical simplex gets a canonical basis vector.
 * Each paired simplex is assigned a vector obtained via XOR of the
 * values of the faces in the boundary (excluding its lower pair).
 *
 * @param ms The Morse sequence context.
 * @param W The Morse sequence.
 */
RefMap::RefMap(MorseSequence& ms, m_sequence W) : MorseFrameBase(ms, W) {
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
            // Critical simplex gets a canonical basis vector
            node_ptr crit_ptr = std::get<node_ptr>(item);
            bitmap b(dim_crit);
            b.set(critToIndex.at(crit_ptr));
            add(crit_ptr, b);
        }
        else if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);

            // Upper simplex gets the zero vector
            add(tau_ptr, bitmap(dim_crit));

            // Compute XOR over boundary of tau, excluding sigma
            node_list bd = ms.boundary(tau_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), sigma_ptr), bd.end());

            bitmap b_sigma(dim_crit);
            for (node_ptr cn : bd) {
                b_sigma ^= get(cn);
            }
            add(sigma_ptr, b_sigma);
        }
        else {
            throw std::runtime_error("Unknown simplex type encountered in Morse sequence.");
        }
    }
}

/**
 * @brief Computes the persistence pairs and essential simplices from the reference map.
 * 
 * The persistence is a list of pairs ((ν, σ), λ) where ν is a critical simplex,
 * σ is the critical simplex paired with ν, and λ is the persistence value (difference in weights).
 * Essential simplices are those with infinite persistence.
 * 
 * @return A pair consisting of:
 *         - A list of essential critical simplices (node_list).
 *         - A vector of persistence pairs (vector of pairs ((node_ptr, node_ptr), int)) sorted by persistence ascending.
 * 
 * @note The persistence result is analogous to the copersistence result computed from the coreference map
 *       when both maps are constructed from the same Morse sequence.
 */

// To test on a copy of the original refmap
// RefMap morse_frame = ref_map;
// morse_frame.persistence();
std::pair<node_list, std::vector<std::pair<node_pair, int>>> RefMap::persistence() {
    node_list essential;
    std::vector<std::pair<node_pair, int>> pairs;

    for (node_ptr sigma : critics) {
        // Compute the boundary vector of sigma
        node_list boundary = ms.boundary(sigma);
        bitmap bsigma(dim_crit);
        for (node_ptr tau : boundary) {
            bsigma ^= get(tau);  // XOR the bitarrays from the boundary
        }

        if (bsigma.any()) {
            // Find the critical simplex nu that maximizes F[tau]
            double max_value = -std::numeric_limits<double>::infinity();
            node_ptr nu;
            node_stack F = ms.get_stack();
            for (size_t i = bsigma.find_first(); i != boost::dynamic_bitset<>::npos; i = bsigma.find_next(i)) {
                node_ptr tau = indexToCrit[i];
                if (F[tau] >= max_value) {
                    max_value = F[tau];
                    nu = tau;
                }
            }

            // Store the persistence pair ((nu, sigma), F[sigma] - F[nu])
            pairs.push_back({{nu, sigma}, F[sigma] - F[nu]});

            // Update all bitarrays to account for the persistence pair
            int index_sigma = critToIndex[sigma];
            int index_nu = critToIndex[nu];
            transform_bitarrays([&](const bitmap& ba) -> bitmap {
                return update_bitarray(ba, index_sigma, index_nu, bsigma);
            });
        }
    }

    // Sort the persistence pairs by increasing persistence
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
