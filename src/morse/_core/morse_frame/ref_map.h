#ifndef REF_MAP_H
#define REF_MAP_H

#include "morse_frame/utils.h"

/**
 * @brief Class computing the reference map from a Morse sequence.
 */
class RefMap : public MorseFrameBase {
public:
    /**
     * @brief Construct a RefMap and compute the corresponding Morse frame.
     *
     * Each critical simplex is assigned a canonical unit vector.
     * Each pair (sigma, tau) is processed by giving tau the zero vector,
     * and sigma the XOR of its boundary faces (excluding sigma itself).
     *
     * @param ms The Morse sequence context (e.g., for boundary/coboundary ops).
     * @param W The Morse sequence (list of critical simplices and paired simplex pairs).
     */
    explicit RefMap(MorseSequence& ms, m_sequence W);

    /**
     * @brief Computes the persistence pairs and essential simplices from the reference map.
     * 
     * @return A pair consisting of:
     *         - A list of essential critical simplices (node_list).
     *         - A vector of persistence pairs (vector of pairs ((node_ptr, node_ptr), int)) sorted by persistence ascending.
     */
    std::pair< node_list, std::vector<std::pair<node_pair, int>> > persistence();
};

#endif  // REF_MAP_H