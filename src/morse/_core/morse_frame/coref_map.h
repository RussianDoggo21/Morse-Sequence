#ifndef COREF_MAP_HPP
#define COREF_MAP_HPP

#include "morse_frame_base.h"

/**
 * @brief Class computing the coreference map from a Morse sequence.
 */

class CorefMap : public MorseFrameBase {
public:
    explicit CorefMap(MorseSequence& ms, m_sequence W);
};

#endif  // COREF_MAP_HPP
