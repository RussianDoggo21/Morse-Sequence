#ifndef UNION_FIND_H
#define UNION_FIND_H

#include "morse_sequence.h" 
#include <boost/dynamic_bitset.hpp>

using bitmap = boost::dynamic_bitset<>;

class UnionFind {
public:
    UnionFind();

    node_ptr _find(const node_ptr& x);
    void _union(const node_ptr& x, const node_ptr& y);
    void add(const node_ptr& x, const bitmap& bitarray_val);

private:
    tsl::robin_map<node_ptr, node_ptr> parent;
    tsl::robin_map<node_ptr, int> size;
    tsl::robin_map<node_ptr, bitmap> bitarray;
    tsl::robin_map<bitmap, node_ptr> bitarray_to_root;
};

#endif
