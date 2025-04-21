#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "simplextree-py/include/simplextree.h"  // Inclusion du complexe simplicial
#include <optional>
#include <variant>

class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  // Constructeur prenant un SimplexTree
    //int dim(const simplex_t& sigma);
    void union_(std::vector<simplex_t>& l, const simplex_t& elmt);
    //void print_simplex(const simplex_t& sigma) const;
    void print_simplex(node_ptr cn);
    //vector<simplex_t> boundary(const simplex_t& sigma);
    vector<node_ptr> boundary(node_ptr cn);
    //vector<simplex_t> coboundary(const simplex_t& sigma);
    vector<node_ptr> coboundary(node_ptr cn);
    //std::optional<simplex_t> find_out(const vector<simplex_t>& B, set<simplex_t>& S);
    std::optional<node_ptr> find_out(vector<node_ptr> B, set<node_ptr> S);
    //vector<simplex_t> tri_dim_decroissant(vector<simplex_t>& L);
    vector<node_ptr> tri_dim_decroissant(vector<node_ptr> L); 
    //vector<simplex_t> tri_dim_croissant(vector<simplex_t>& L);
    vector<node_ptr> tri_dim_croissant(vector<node_ptr> L);
    //vector<simplex_t> simplices(std::optional<int> p = std::nullopt) const;
    vector<node_ptr> simplices(std::optional<int> p = std::nullopt) const;
    //std::pair<std::vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>>, int> morse_seq_dec(const SimplexTree& st);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> morse_seq_dec(const SimplexTree& st);
    //std::pair<std::vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>>, int> morse_seq_crois(const SimplexTree& K_init);
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> morse_seq_crois(const SimplexTree& st);

private:
    const SimplexTree& simplex_tree;  // Référence au complexe simplicial
};

#endif // MORSE_SEQUENCE_H
