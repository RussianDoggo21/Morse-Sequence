#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "../../../simplextree-py/include/simplextree.h"  // Inclusion of the library SimplexTree
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

#include <tsl/robin_map.h>
#include <boost/dynamic_bitset.hpp>

using node_pair = std::pair<node_ptr, node_ptr>;
using m_sequence = std::vector<std::variant<node_ptr, node_pair>>;
using node_list = std::vector<node_ptr>;
using bitmap = boost::dynamic_bitset<>;
using m_frame = tsl::robin_map<node_ptr, bitmap>;
using m_frame0 = tsl::robin_map<node_ptr, node_list>;
using node_index_map = tsl::robin_map<node_ptr, std::size_t>;
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

    // Morse Frames (bitmap implementation)
    node_index_map generate_critical_index_map(const m_sequence& W);
    m_frame reference_map(const m_sequence& W, const node_index_map& critical_index_map);
    m_frame coreference_map(const m_sequence& W, const node_index_map& critical_index_map);
    void print_m_frame(const m_frame& map, const m_sequence& W, const node_index_map& critical_index_map);

    // Morse Frames (simplex implementation)
    m_frame0 reference_map0(const m_sequence& W);
    m_frame0 coreference_map0(const m_sequence& W);
    void print_m_frame0(m_frame0& map, const m_sequence& W);
    

private:
    const SimplexTree& simplex_tree;  // Reference to the simplicial complex given in input
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, const node_ptr& s_ptr);
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const tsl::robin_map<node_ptr, int>& F);
    node_list sym_diff(const node_list& A, const node_list& B);
    void print_bitmap(const bitmap& bm, const m_sequence& W, const node_index_map& critical_index_map) const;
};

#endif 