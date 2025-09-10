#ifndef MORSE_FRAME_BASE_H
#define MORSE_FRAME_BASE_H

#include "morse_frame/union_find_mf.h"

#include <limits>

using m_frame = tsl::robin_map<node_ptr, bitmap>;
using node_index_map = tsl::robin_map<node_ptr, std::size_t>;
using index_node_map = tsl::robin_map<std::size_t, node_ptr>;


struct FullBitArrayElement {
    bool is_critical;  
    tsl::robin_map<node_ptr, bitmap> critical_simplex;  
    std::pair<std::pair<node_ptr, bitmap>, std::pair<node_ptr, bitmap>> pair;
};
typedef std::vector<FullBitArrayElement> FullBitArray;


/**
 * @brief Base class for computing Morse frames.
 * 
 * Provides utility functions to manage critical simplices and
 * manipulate bitmaps associated to each node.
 */
class Utils : public UnionFindMF {
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

    /**
     * @brief Update the bitarray during persistence or copersistence computation.
     *
     * @param ba The bitarray to update (passed by reference).
     * @param indexsigma Index of the critical simplex sigma in the bitarray.
     * @param indexnu Index of the critical simplex nu in the bitarray.
     * @param bsigma The bitarray representing the boundary of sigma.
     *
     * @return The updated bitarray.
     */
    bitmap update_bitarray(const bitmap& ba, size_t indexsigma, size_t indexnu, const bitmap& bsigma) const;

    FullBitArray full_bitarray;

public:
    /**
     * @brief Constructor initializes the Morse frame base.
     * 
     * @param ms MorseSequence instance reference.
     * @param W Sequence of critical simplices and pairs.
     */
    Utils(MorseSequence& ms, const m_sequence& W);

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

    void print_persistence_results(const std::pair<node_list, std::vector<std::pair<node_pair, int>>>& results) const;

    FullBitArray get_full_bitarray() const;

    void print_full_bitarray() const;
};

#endif // MORSE_FRAME_BASE_H
