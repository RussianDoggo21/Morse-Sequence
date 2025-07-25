#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "../../../../simplextree-py/include/simplextree.h"  // Inclusion of the library SimplexTree
#include <optional>
#include <variant>
#include <list>
#include <vector>
#include <unordered_map>
#include <algorithm> 
#include <tuple>
#include <cstddef>
#include <functional>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

#include <tsl/robin_map.h>

using node_pair = std::pair<node_ptr, node_ptr>;
using m_sequence = std::vector<std::variant<node_ptr, node_pair>>;
using node_list = std::vector<node_ptr>;
using SimplexList = std::vector<simplex_t>;  
using simplex_t = SimplexTree::simplex_t;


class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  // Constructor

    // Auxiliary functions
    const SimplexTree& get_simplex_tree();
    node_list boundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);
    node_list boundary(const node_ptr& cn); 
    node_list coboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);
    node_list coboundary(const node_ptr& cn);
    int nbboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);
    int nbboundary(const node_ptr& cn);
    int nbcoboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);
    int nbcoboundary(const node_ptr& cn);
    node_list simplices(std::optional<int> p) const;

    // Morse Sequences - F-Sequences
    std::pair<m_sequence, int> increasing(const SimplexTree& st);
    std::pair<m_sequence, int> decreasing(const SimplexTree& st);
    std::pair<m_sequence, int> Max(const node_list& S, const tsl::robin_map<node_ptr, int>& F);
    std::pair<m_sequence, int> Min(const node_list& S, const tsl::robin_map<node_ptr, int>& F);
    void print_morse_sequence(const std::pair<m_sequence, int>& result, bool n_crit = false);
    

private:
    const SimplexTree& simplex_tree;  // Reference to the simplicial complex given in input
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, const node_ptr& s_ptr);
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const tsl::robin_map<node_ptr, int>& F);
};

#endif 