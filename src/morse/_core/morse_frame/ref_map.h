#ifndef REF_MAP_HPP
#define REF_MAP_HPP

#include "morse_frame_base.h"

/**
 * @brief Class computing the reference map from a Morse sequence.
 */

class RefMap : public MorseFrameBase {
public:
    explicit RefMap(MorseSequence& ms, m_sequence W); 
};

#endif  // REF_MAP_HPP
