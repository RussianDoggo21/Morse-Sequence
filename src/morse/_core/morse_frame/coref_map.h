#ifndef COREF_MAP_H
#define COREF_MAP_H

#include "morse_frame/morse_frame_base.h"

/**
 * @brief Class computing the *coreference map* from a Morse sequence.
 *
 * In contrast to RefMap (reference map), the computation here proceeds in reverse order
 * and uses *coboundaries* rather than boundaries.
 */
class CorefMap : public MorseFrameBase {
public:
    /**
     * @brief Construct a CorefMap and compute the corresponding Morse frame.
     *
     * The Morse sequence is processed in reverse.
     * Critical simplices get unit vectors, while paired ones follow
     * an XOR logic over their coboundaries.
     *
     * @param ms The Morse sequence context (e.g., for boundary/coboundary ops).
     * @param W The Morse sequence (list of critical simplices and paired simplex pairs).
     */
    explicit CorefMap(MorseSequence& ms, m_sequence W);
};

#endif  // COREF_MAP_H