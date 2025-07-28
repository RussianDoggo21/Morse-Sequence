#ifndef REF_MAP_H
#define REF_MAP_H

#include "morse_frame_base.h"

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

    std::pair<node_list, std::vector<node_pair>> persistence();
};

#endif  // REF_MAP_H