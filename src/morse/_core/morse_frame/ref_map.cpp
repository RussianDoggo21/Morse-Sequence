#include "ref_map.h"

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
