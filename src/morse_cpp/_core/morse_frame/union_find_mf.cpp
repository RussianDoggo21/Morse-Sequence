#include "morse_frame/union_find_mf.h"

UnionFindMF::UnionFindMF() {}

/**
 * @brief Find the root of a given node with path compression.
 * 
 * This function compresses the path by making each visited node
 * point directly to the root.
 * 
 * @param x The node_ptr to find the root of.
 * @return The root node_ptr of the connected component.
 */
node_ptr UnionFindMF::_find(const node_ptr& x) {
    if (parent[x] != x) {
        parent[x] = _find(parent[x]);  // Path compression
    }
    return parent[x];
}

/**
 * @brief Merge the components of x and y using union by size.
 * @param x First node_ptr.
 * @param y Second node_ptr.
 */
void UnionFindMF::_union(const node_ptr& x, const node_ptr& y) {
    node_ptr x_root = _find(x);
    node_ptr y_root = _find(y);
    if (x_root == y_root) return;  // Already in the same set

    // Ensure x_root is the larger set
    if (size[x_root] < size[y_root]) std::swap(x_root, y_root);

    parent[y_root] = x_root;
    size[x_root] += size[y_root];
}

/**
 * @brief Add a simplex with its bitmap.
 * 
 * If the bitmap already exists, merge the simplex into the corresponding equivalence class.
 * Otherwise, create a new class for it.
 * 
 * @param x The simplex (node_ptr) to add.
 * @param bitarray_val The bitmap associated to this simplex.
 */
void UnionFindMF::add(const node_ptr& x, const bitmap& bitarray_val) {
    auto it = bitarray_to_root.find(bitarray_val);

    // Case 1: An equivalence class for bitarray_val already exists
    if (it != bitarray_to_root.end()) {
        // Retrieve the representative (root) of the class of bitarray_val
        node_ptr existing_root = it->second;

        // Subcase A: 
        // if parent.count(x) != 0, it means that x is the root of an existing class 
        // different from bitarray_val
        // we merge the existing class of x with the new class represented by bitarray_val
        if (parent.count(x)) {
            _union(x, existing_root);
        }

        // Subcase B: x is new — attach it directly to the existing class
        else {
            parent[x] = existing_root;
            size[existing_root]++;
        }
        return;
    }

    // Case 2: No equivalence class exists for bitarray_val
    // Create a new component with x as its own representative
    parent[x] = x;
    size[x] = 1;
    bitarray[x] = bitarray_val;
    bitarray_to_root[bitarray_val] = x;
}

/**
 * @brief Get the bitmap associated to the representative of x.
 * 
 * Throws an exception if the bitmap is not found.
 * 
 * @param x The simplex (node_ptr) to query.
 * @return The bitmap associated with the root of x.
 */
bitmap UnionFindMF::get(const node_ptr& x) {
    auto root = this->_find(x);
    auto it = bitarray.find(root);
    if (it != bitarray.end()) {
        return it->second;
    }
    throw std::runtime_error("Simplex not found.");
}

/**
 * @brief Return the map linking root nodes to their bitmap.
 * @return Map from node_ptr to bitmap.
 */
tsl::robin_map<node_ptr, bitmap> UnionFindMF::get_bitarray() const {
    return this->bitarray;
}

/**
 * @brief Applies a transformation to all bitarrays stored in the structure.
 *
 * This function iterates over all bitarrays associated with Union-Find roots,
 * applies a user-defined transformation, and performs component merging
 * if the transformed bitarray already exists in the structure.
 *
 * If a transformed bitarray is already present, the corresponding components
 * are merged via `_union`. If not, the current component is updated with the new
 * bitarray.
 *
 * @param func A function taking a `const bitmap&` and returning a `bitmap`,
 *             representing the transformation to apply.
 */
void UnionFindMF::transform_bitarrays(const std::function<bitmap(const bitmap&)>& func) {
    // Snapshot of the original bitarrays to avoid modifying the map during iteration
    std::vector<bitmap> originals;
    originals.reserve(bitarray_to_root.size());
    for (const auto& [ba, _] : bitarray_to_root) {
        originals.push_back(ba);
    }

    for (const bitmap& old_ba : originals) {
        auto old_root_it = bitarray_to_root.find(old_ba);
        if (old_root_it == bitarray_to_root.end()) continue;
        node_ptr old_root = old_root_it->second;

        bitmap new_ba = func(old_ba);

        if (new_ba == old_ba) continue;

        auto new_root_it = bitarray_to_root.find(new_ba);
        if (new_root_it != bitarray_to_root.end()) {
            node_ptr new_root = new_root_it->second;
            if (old_root != new_root) {
                _union(old_root, new_root);
                node_ptr merged_root = _find(old_root);

                bitarray[merged_root] = new_ba;
                bitarray_to_root[new_ba] = merged_root;

                if (old_root != merged_root) {
                    bitarray.erase(old_root);
                }
                if (new_root != merged_root) {
                    bitarray.erase(new_root);
                }
            }
            bitarray_to_root.erase(old_ba);
        } else {
            bitarray[old_root] = new_ba;
            bitarray_to_root[new_ba] = old_root;
            bitarray_to_root.erase(old_ba);
        }
    }
}
