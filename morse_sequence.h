#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "simplextree-py/include/simplextree.h"  // Inclusion du complexe simplicial
#include <optional>
#include <variant>

class MorseSequence {
public:
    explicit MorseSequence(const SimplexTree& st);  // Constructeur prenant un SimplexTree
    vector<node_ptr> boundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    vector<node_ptr> coboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    int nbboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    int nbcoboundary(node_ptr cn, const unordered_map<node_ptr, bool>& S);
    vector<node_ptr> simplices(std::optional<int> p) const;
    std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> Max(const vector<node_ptr>& S, const unordered_map<simplex_t, int>& F);

private:
    const SimplexTree& simplex_tree;  // Référence au complexe simplicial
};

#endif // MORSE_SEQUENCE_H
