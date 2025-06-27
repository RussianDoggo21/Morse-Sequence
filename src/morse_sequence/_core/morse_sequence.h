#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "../../../simplextree-py/include/simplextree.h"  // Inclusion of the library SimplexTree
#include <optional>
#include <variant>
#include <list>
#include <tsl/robin_map.h>

using node_pair = std::pair<node_ptr, node_ptr>;
using m_sequence = std::vector<std::variant<node_ptr, node_pair>>;
using node_list = std::vector<node_ptr>;
using m_frame = std::unordered_map<node_ptr, node_list>;
using node_map = std::unordered_map<node_ptr,node_ptr>;
using simplex_t = SimplexTree::simplex_t;
/*
typedef std::size_t idx_t;
using simplex_t = vector< idx_t >; 
*/

class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  
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
    //node_list get_node_list(const std::list<simplex_t>& py_list) const;
    std::pair<m_sequence, int> increasing(const SimplexTree& st);
    std::pair<m_sequence, int> decreasing(const SimplexTree& st);
    std::pair<m_sequence, int> Max(const node_list& S, const unordered_map<node_ptr, int>& F);
    std::pair<m_sequence, int> Min(const node_list& S, const std::unordered_map<node_ptr, int>& F);
    void print_morse_sequence(const std::pair<m_sequence, int>& result, bool n_crit = false);
    m_frame reference_map(const m_sequence& W);
    m_frame coreference_map(const m_sequence& W);
    void print_m_frame(m_frame& map, const m_sequence& W);

private:
    const SimplexTree& simplex_tree;  // Reference to the simplicial complex given in input
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, const node_ptr& s_ptr);
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T,const node_list& simplex_list, node_ptr s_ptr, const std::unordered_map<node_ptr, int>& F);
    node_list sym_diff(const node_list& A, const node_list& B);
    enum class GammaMode { Reference, Coreference };
    node_list Gamma(const node_ptr& sigma_ptr, const node_map& sigma2tau, const node_map& tau2sigma, m_frame& cache, GammaMode mode);
};

#endif 
