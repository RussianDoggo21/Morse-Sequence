#include "coref_map.h"

CorefMap::CorefMap(MorseSequence& ms, m_sequence W) : MorseFrameBase(ms, W) {
    // Compute the Morse frame for the coreference map (in reverse order)
    for (auto it = W.rbegin(); it != W.rend(); ++it) {
        const auto& item = *it;

        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit_ptr = std::get<node_ptr>(item);
            bitmap b(dim_crit);
            b.set(critToIndex.at(crit_ptr));
            add(crit_ptr, b);
        }
        else if (std::holds_alternative<node_pair>(item)) {
            auto [sigma_ptr, tau_ptr] = std::get<node_pair>(item);

            add(sigma_ptr, bitmap(dim_crit));  // σ gets zero vector

            node_list cbd = ms.coboundary(sigma_ptr);
            cbd.erase(std::remove(cbd.begin(), cbd.end(), tau_ptr), cbd.end());

            bitmap b_tau(dim_crit);
            for (node_ptr cn : cbd) {
                b_tau ^= get(cn);  // XOR with existing values
            }
            add(tau_ptr, b_tau);
        }
        else {
            throw std::runtime_error("Unknown simplex type encountered in Morse sequence.");
        }
    }
}
