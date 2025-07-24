#ifndef MORSE_FRAME_H
#define MORSE_FRAME_H

#include "union_find.h"

using m_frame = tsl::robin_map<node_ptr, bitmap>;
using m_frame0 = tsl::robin_map<node_ptr, node_list>;
using node_index_map = tsl::robin_map<node_ptr, std::size_t>;
using index_node_map = tsl::robin_map<std::size_t, node_ptr>;

class RefMap : public UnionFind{
    public:
        explicit RefMap(MorseSequence& ms);

        // Morse Frames (bitmap implementation)
        m_frame reference_map(const m_sequence& W);
        m_frame coreference_map(const m_sequence& W);
        void print_m_frame(const m_frame& map, const m_sequence& W);

        // Morse Frames (simplex implementation)
        m_frame0 reference_map0(const m_sequence& W);
        m_frame0 coreference_map0(const m_sequence& W);
        void print_m_frame0(m_frame0& map, const m_sequence& W);
    
    private:

        // Attributes of the class
        MorseSequence& ms; // Reference to the Morse Sequence
        const SimplexTree& simplex_tree; // Reference to the simplex tree
        node_list critics; // Critical simplices from ms
        node_index_map critToIndex; // Map critical -> index
        index_node_map indexToCrit; // Map index -> critical
        int dim_crit; // Dimension of the bitarray (number of critical simplices)

        // Auxiliary functions
        void generate_attributes(const m_sequence& W);
        node_list sym_diff(const node_list& A, const node_list& B);
        void print_bitmap(const bitmap& bm, const m_sequence& W) const;
};

#endif