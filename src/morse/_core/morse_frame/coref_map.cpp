#include "coref_map.h"

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
CorefMap::CorefMap(MorseSequence& ms, m_sequence W) : MorseFrameBase(ms, W) {
    for (auto it = W.rbegin(); it != W.rend(); ++it) {
        const auto& item = *it;

        if (std::holds_alternative<node_ptr>(item)) {
            // Critical simplex: assign canonical unit vector
            node_ptr crit_ptr = std::get<node_ptr>(item);
            bitmap b(dim_crit);
            b.set(critToIndex.at(crit_ptr));
            add(crit_ptr, b);
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
        }

        else {
            throw std::runtime_error("Unknown simplex type encountered in Morse sequence.");
        }
    }
}
