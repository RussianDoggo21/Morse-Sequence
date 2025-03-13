#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "simplextree-py/include/simplextree.h"  // Inclusion du complexe simplicial
#include <optional>
#include <variant>

class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  // Constructeur prenant un SimplexTree
    int dim(const simplex_t& sigma);
    void union_(std::vector<simplex_t>& l, const simplex_t& elmt);
    void print_simplex(const simplex_t& sigma) const;
    vector<simplex_t> boundary(const simplex_t& sigma);
    vector<simplex_t> cofaces(const simplex_t& sigma);
    vector<simplex_t> coboundary(const simplex_t& sigma);
    std::optional<simplex_t> find_out(const vector<simplex_t>& coboundary, vector<simplex_t>& S);
    vector<simplex_t> tri_dim_decroissant(vector<simplex_t>& L);
    vector<simplex_t> tri_dim_croissant(vector<simplex_t>& L);
    vector<simplex_t> simplices(std::optional<int> p = std::nullopt) const;
    std::pair<std::vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>>, int> morse_seq_dec(const SimplexTree& st);
    std::pair<std::vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>>, int> morse_seq_crois(const SimplexTree& K_init);

private:
    const SimplexTree& simplex_tree;  // Référence au complexe simplicial
};

#endif // MORSE_SEQUENCE_H
