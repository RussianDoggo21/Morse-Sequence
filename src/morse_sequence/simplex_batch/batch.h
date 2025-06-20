#ifndef SIMPLEX_BATCH_H
#define SIMPLEX_BATCH_H

#include "../../../simplextree-py/include/simplextree.h"
#include <vector>

// To reduce compilation time
// Declaration of pybind11::object in batch.cpp
namespace pybind11 { class object; } 

// Structure to reduce the conversion between python objects and C++ objects
struct SimplexBatch{

    std::vector<simplex_t> simplices;   // Python faces (original order)
    std::vector<node_ptr>  nodes;       // aligned node_ptrs
    std::vector<int>       weights;     // aligned weights (empty for S)

    static SimplexBatch from_python(const pybind11::object& obj, SimplexTree& tree);
};

// Generates the unordered_map<node_ptr, int> needed when calling the Min and Max function of MorseSequence  
inline std::unordered_map<node_ptr,int> batch_to_weight_map(const SimplexBatch& batch){
    std::unordered_map<node_ptr,int> out;
    out.reserve(batch.nodes.size());
    for (std::size_t i = 0; i < batch.nodes.size(); ++i)
        out.emplace(batch.nodes[i],
                    batch.weights.empty() ? 0 : batch.weights[i]);
    return out;
}

#endif /* SIMPLEX_BATCH_H */
