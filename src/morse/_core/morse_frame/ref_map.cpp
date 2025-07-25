#include "ref_map.h"

RefMap::RefMap(MorseSequence& ms, m_sequence W) : MorseFrameBase(ms, W){        
        
    // Compute the Morse frame
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit_ptr = std::get<node_ptr>(item);
            bitmap b(dim_crit);
            b.set(critToIndex.at(crit_ptr));
            add(crit_ptr, b);
        }

        else if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);

            add(tau_ptr, bitmap(dim_crit));  // τ gets zero vector

            node_list bd = ms.boundary(tau_ptr);
            bd.erase(std::remove(bd.begin(), bd.end(), sigma_ptr), bd.end());

            bitmap b_sigma(dim_crit);
            for (node_ptr cn : bd) {
                b_sigma ^= get(cn);  // XOR with existing values
            }
            add(sigma_ptr, b_sigma);
        }

        else {
            throw std::runtime_error("Unknown simplex type encountered in Morse sequence.");
        }
    }
}