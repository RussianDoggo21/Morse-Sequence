#ifndef MORSE_FRAME_BASE_H
#define MORSE_FRAME_BASE_H

#include "union_find.h"

using m_frame = tsl::robin_map<node_ptr, bitmap>;
using node_index_map = tsl::robin_map<node_ptr, std::size_t>;
using index_node_map = tsl::robin_map<std::size_t, node_ptr>;

class MorseFrameBase : public UnionFindMF {
protected:
    MorseSequence& ms;
    const SimplexTree& simplex_tree;
    m_sequence W;
    std::vector<node_ptr> critics;
    node_index_map critToIndex;
    std::unordered_map<size_t, node_ptr> indexToCrit;
    size_t dim_crit;

public:
    MorseFrameBase(MorseSequence& ms, const m_sequence& W);

    void print_bitmap(const bitmap& bm, const m_sequence& W) const;
    void print_m_frame(const m_sequence& W) const;
};

#endif // MORSE_MAP_BASE_H
