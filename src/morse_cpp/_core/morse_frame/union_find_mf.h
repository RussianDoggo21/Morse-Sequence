#ifndef UNION_FIND_H
#define UNION_FIND_H

#include "morse_sequence/morse_sequence.h" 
#include <boost/dynamic_bitset.hpp>

using bitmap = boost::dynamic_bitset<>;

/**
 * @brief Union-Find structure for Morse Frame computation.
 * 
 * Maintains connected components of simplices (node_ptr)
 * and associates each root with a bitmap representing a 
 * Morse frame vector.
 */
class UnionFindMF {
public:
    UnionFindMF();

    /**
     * @brief Find the representative (root) of a simplex with path compression.
     * @param x The simplex (node_ptr) to find the root of.
     * @return The root node_ptr representing the connected component.
     */
    node_ptr _find(const node_ptr& x);


    /**
     * @brief Merge the connected components of two simplices.
     * @param x First simplex.
     * @param y Second simplex.
     */
    void _union(const node_ptr& x, const node_ptr& y);


    /**
     * @brief Add a simplex with a given bitmap.
     * 
     * If the same bitmap already exists, the simplex is merged
     * with the class of the existing one. Otherwise, it becomes
     * the root of a new component.
     * 
     * @param x The simplex to add.
     * @param bitarray_val The bitmap associated to this simplex.
     */
    void add(const node_ptr& x, const bitmap& bitarray_val);


    /**
     * @brief Return the bitmap associated with the component of x.
     * @param x The simplex to query.
     * @return The bitmap associated to the root of x.
     */
    bitmap get(const node_ptr& x);


    /**
     * @brief Return a copy of the internal map from root nodes to bitmaps.
     * @return A map from node_ptr to bitmap.
     */
    tsl::robin_map<node_ptr, bitmap> get_bitarray() const;

    /**
     * @brief Applies a transformation to all bitarrays stored in the structure.
     * @param func A function taking a `const bitmap&` and returning a `bitmap`,
     *             representing the transformation to apply.
     */
    void transform_bitarrays(const std::function<bitmap(const bitmap&)>& func);
    

private:
    // Disjoint set structure: maps each node to its parent
    tsl::robin_map<node_ptr, node_ptr> parent;

    // Stores the size of each set (for union by size)
    tsl::robin_map<node_ptr, int> size;

    // Maps root node to its associated bitmap
    tsl::robin_map<node_ptr, bitmap> bitarray;

    // Maps bitmaps to their associated root (used to detect duplicate bitmaps)
    tsl::robin_map<bitmap, node_ptr> bitarray_to_root;
};

#endif
