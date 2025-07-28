#ifndef MORSE_FRAME_BASE_H
#define MORSE_FRAME_BASE_H

#include "union_find_mf.h"

#include <limits>

using m_frame = tsl::robin_map<node_ptr, bitmap>;
using node_index_map = tsl::robin_map<node_ptr, std::size_t>;
using index_node_map = tsl::robin_map<std::size_t, node_ptr>;

/**
 * @brief Base class for computing Morse frames.
 * 
 * Provides utility functions to manage critical simplices and
 * manipulate bitmaps associated to each node.
 */
class MorseFrameBase : public UnionFindMF {
protected:
    MorseSequence& ms;
    const SimplexTree& simplex_tree;
    m_sequence W;

    // Critical simplices
    std::vector<node_ptr> critics;

    // Mapping from critical simplex to index in bitmap
    node_index_map critToIndex;
    // Mapping from index in bitmap to critical simplex
    tsl::robin_map<size_t, node_ptr> indexToCrit;

    // Number of critical simplices
    size_t dim_crit;

public:
    /**
     * @brief Constructor initializes the Morse frame base.
     * 
     * @param ms MorseSequence instance reference.
     * @param W Sequence of critical simplices and pairs.
     */
    MorseFrameBase(MorseSequence& ms, const m_sequence& W);

    /**
     * @brief Print the bitmap as a list of critical simplices.
     * @param bm The bitmap to print.
     * @param W The sequence of critical simplices and pairs.
     */
    void print_bitmap(const bitmap& bm, const m_sequence& W) const;

    /**
     * @brief Print all entries in the Morse frame with their associated bitmaps.
     * @param W The sequence of critical simplices and pairs.
     */
    void print_m_frame(const m_sequence& W);
};

#endif // MORSE_FRAME_BASE_H
