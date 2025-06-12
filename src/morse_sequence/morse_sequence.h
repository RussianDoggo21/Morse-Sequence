#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "../../simplextree-py/include/simplextree.h"  // Inclusion of the library SimplexTree
#include <optional>
#include <variant>

class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  
    const SimplexTree& get_simplex_tree();
    vector<node_ptr> boundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    vector<node_ptr> boundary(node_ptr cn); 
    vector<node_ptr> coboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    vector<node_ptr> coboundary(node_ptr cn);
    int nbboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    int nbboundary(node_ptr cn);
    int nbcoboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    int nbcoboundary(node_ptr cn);
    vector<node_ptr> simplices(std::optional<int> p) const;
    node_ptr find_out(const std::unordered_map<node_ptr, bool>& T, const std::vector<node_ptr>& simplex_list, std::string order, node_ptr s_ptr);
    node_ptr find_out(const std::unordered_map<node_ptr, bool>& T,const std::vector<node_ptr>& simplex_list, node_ptr s_ptr, const std::unordered_map<node_ptr, int>& F);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> increasing(const SimplexTree& st);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> decreasing(const SimplexTree& st);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> Max(const vector<node_ptr>& S, const unordered_map<node_ptr, int>& F);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> Min(const std::vector<node_ptr>& S, const std::unordered_map<node_ptr, int>& F);

private:
    const SimplexTree& simplex_tree;  // Reference to the simplicial complex given in input
};

#endif 
