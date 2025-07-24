#include "union_find.h"

UnionFind::UnionFind() {}

node_ptr UnionFind::_find(const node_ptr& x) {
    if (parent[x] != x) {
        parent[x] = _find(parent[x]);  // Path compression
    }
    return parent[x];
}

void UnionFind::_union(const node_ptr& x, const node_ptr& y) {
    node_ptr x_root = _find(x);
    node_ptr y_root = _find(y);
    if (x_root == y_root) return;

    if (size[x_root] < size[y_root]) std::swap(x_root, y_root);

    parent[y_root] = x_root;
    size[x_root] += size[y_root];
}

void UnionFind::add(const node_ptr& x, const bitmap& bitarray_val) {
    auto it = bitarray_to_root.find(bitarray_val);
    if (it != bitarray_to_root.end()) {
        node_ptr existing_root = it->second;

        if (parent.count(x)) {
            _union(x, existing_root);
        } else {
            parent[x] = existing_root;
            size[existing_root]++;
        }
        return;
    }

    parent[x] = x;
    size[x] = 1;
    bitarray[x] = bitarray_val;
    bitarray_to_root[bitarray_val] = x;
}
