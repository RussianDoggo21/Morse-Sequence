#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "../../simplextree-py/include/simplextree.h"  // Inclusion of the library SimplexTree
#include <optional>
#include <variant>
using node_pair = std::pair<node_ptr, node_ptr>;
using m_sequence = std::vector<std::variant<node_ptr, node_pair>>;
using node_list = std::vector<node_ptr>;
using morse_frame = std::unordered_map<node_ptr, node_list>;
using node_map = std::unordered_map<node_ptr,node_ptr>;


class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  
    const SimplexTree& get_simplex_tree();
    node_list boundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    node_list boundary(node_ptr cn); 
    node_list coboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    node_list coboundary(node_ptr cn);
    int nbboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    int nbboundary(node_ptr cn);
    int nbcoboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    int nbcoboundary(node_ptr cn);
    node_list simplices(std::optional<int> p) const;
    std::pair<m_sequence, int> increasing(const SimplexTree& st);
    std::pair<m_sequence, int> decreasing(const SimplexTree& st);
    std::pair<m_sequence, int> Max(const node_list& S, const unordered_map<node_ptr, int>& F);
    std::pair<m_sequence, int> Min(const node_list& S, const std::unordered_map<node_ptr, int>& F);
    void print_morse_sequence(std::pair<m_sequence, int> result, bool n_crit = false);
    morse_frame reference_map(const m_sequence& W);
    morse_frame coreference_map(const m_sequence& W);
    void print_morse_frame(morse_frame& map, const m_sequence& W);

private:
    const SimplexTree& simplex_tree;  // Reference to the simplicial complex given in input
    node_ptr find_out(const std::unordered_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, node_ptr s_ptr);
    node_ptr find_out(const std::unordered_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const std::unordered_map<node_ptr, int>& F);
    node_list sym_diff(const node_list& A, const node_list& B);
    enum class GammaMode { Reference, Coreference };
    node_list Gamma(node_ptr sigma_ptr, const node_map& sigma2tau, const node_map& tau2sigma, morse_frame& cache, GammaMode mode);
};

#endif 
