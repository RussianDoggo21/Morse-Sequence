#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "simplextree-py/include/simplextree.h"  // Inclusion du complexe simplicial
#include <optional>
#include <variant>

class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  // Constructeur prenant un SimplexTree
    void union_(std::vector<simplex_t>& l, const simplex_t& elmt);
    void print_simplex(node_ptr cn);
    vector<node_ptr> boundary(node_ptr cn);
    vector<node_ptr> coboundary(node_ptr cn);
    std::optional<node_ptr> find_out(vector<node_ptr> B, set<node_ptr> S);
    vector<node_ptr> tri_dim_decroissant(vector<node_ptr> L); 
    vector<node_ptr> tri_dim_croissant(vector<node_ptr> L);
    vector<node_ptr> simplices(std::optional<int> p = std::nullopt) const;
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> morse_seq_dec(const SimplexTree& st);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> morse_seq_crois(const SimplexTree& st);

private:
    const SimplexTree& simplex_tree;  // Référence au complexe simplicial
};

#endif // MORSE_SEQUENCE_H
