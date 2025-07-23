#ifndef MORSE_FRAME_H
#define MORSE_FRAME_H

#include "morse_sequence.h"

class MorseFrame{
    public:
        explicit MorseFrame(MorseSequence& ms);

        // Morse Frames (bitmap implementation)
        node_index_map generate_critical_index_map(const m_sequence& W);
        m_frame reference_map(const m_sequence& W, const node_index_map& critical_index_map);
        m_frame coreference_map(const m_sequence& W, const node_index_map& critical_index_map);
        void print_m_frame(const m_frame& map, const m_sequence& W, const node_index_map& critical_index_map);

        // Morse Frames (simplex implementation)
        m_frame0 reference_map0(const m_sequence& W);
        m_frame0 coreference_map0(const m_sequence& W);
        void print_m_frame0(m_frame0& map, const m_sequence& W);
    
    private:
        MorseSequence& ms; // Reference to the Morse Sequence
        const SimplexTree& simplex_tree; // Reference to the simplex tree
        node_list sym_diff(const node_list& A, const node_list& B);
        void print_bitmap(const bitmap& bm, const m_sequence& W, const node_index_map& critical_index_map) const;
};

#endif