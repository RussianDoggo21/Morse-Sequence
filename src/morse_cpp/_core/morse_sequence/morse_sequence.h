#ifndef MORSE_SEQUENCE_H
#define MORSE_SEQUENCE_H

#include "simplextree.h"  // Inclusion of the library SimplexTree
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
using node_stack = tsl::robin_map<node_ptr, int>;

/**
 * @class MorseSequence
 * @brief Class for generating and manipulating Morse sequences on a simplicial complex.
 */
class MorseSequence {
public:
    /**
     * @brief Constructor.
     * @param st Reference to a SimplexTree object representing the simplicial complex.
     */
    explicit MorseSequence(const SimplexTree& st);


    /**
     * @brief Get the underlying simplex tree.
     * @return Reference to the SimplexTree.
     */
    const SimplexTree& get_simplex_tree();


    /**
     * @brief Compute the boundary of a simplex using a filtration.
     * @param cn Pointer to the simplex.
     * @param S A map indicating whether each simplex is included in the filtration.
     * @return List of boundary simplexes.
     */
    node_list boundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);


    /**
     * @brief Compute the boundary of a simplex.
     * @param cn Pointer to the simplex.
     * @return List of boundary simplexes.
     */
    node_list boundary(const node_ptr& cn);


    /**
     * @brief Compute the coboundary of a simplex using a filtration.
     * @param cn Pointer to the simplex.
     * @param S A map indicating whether each simplex is included in the filtration.
     * @return List of coboundary simplexes.
     */
    node_list coboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);


    /**
     * @brief Compute the coboundary of a simplex.
     * @param cn Pointer to the simplex.
     * @return List of coboundary simplexes.
     */
    node_list coboundary(const node_ptr& cn);
    

    /**
     * @brief Count the number of boundary faces using a filtration.
     * @param cn Pointer to the simplex.
     * @param S A map indicating whether each simplex is included in the filtration.
     * @return Number of boundary faces.
     */
    int nbboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);


    /**
     * @brief Count the number of boundary faces.
     * @param cn Pointer to the simplex.
     * @return Number of boundary faces.
     */
    int nbboundary(const node_ptr& cn);


    /**
     * @brief Count the number of coboundary cofaces using a filtration.
     * @param cn Pointer to the simplex.
     * @param S A map indicating whether each simplex is included in the filtration.
     * @return Number of coboundary cofaces.
     */
    int nbcoboundary(const node_ptr& cn, const tsl::robin_map<node_ptr, bool>& S);


    /**
     * @brief Count the number of coboundary cofaces.
     * @param cn Pointer to the simplex.
     * @return Number of coboundary cofaces.
     */
    int nbcoboundary(const node_ptr& cn);


    /**
     * @brief Get the simplices of the simplicial complex.
     * @param p Optional dimension to filter simplices.
     * @return List of simplices.
     */
    node_list simplices(std::optional<int> p) const;


    /**
     * @brief Build an increasing Morse sequence.
     * @param st The input simplicial complex.
     * @return Pair of Morse sequence and number of critical simplices.
     */
    std::pair<m_sequence, int> increasing();


    /**
     * @brief Build a decreasing Morse sequence.
     * @param st The input simplicial complex.
     * @return Pair of Morse sequence and number of critical simplices.
     */
    std::pair<m_sequence, int> decreasing();


    /**
     * @brief Compute a maximal F-sequence based on a function F.
     * @param S List of simplices.
     * @return Pair of Morse sequence and number of critical simplices.
     */
    //std::pair<m_sequence, int> Max(const node_list& S, const node_stack& F);
    std::pair<m_sequence, int> Max(const node_list& S);


    /**
     * @brief Compute a minimal F-sequence based on a function F.
     * @param S List of simplices.
     * @return Pair of Morse sequence and number of critical simplices.
     */
    //std::pair<m_sequence, int> Min(const node_list& S, const node_stack& F);
    std::pair<m_sequence, int> Min(const node_list& S);


    /**
     * @brief Print the Morse sequence and optionally the number of critical simplices.
     * @param result Pair of Morse sequence and number of critical simplices.
     * @param n_crit Whether to print the number of critical simplices.
     */
    void print_morse_sequence(const std::pair<m_sequence, int>& result, bool n_crit = false);

    /**
     * @brief Get the node_stack values associated to the simplicial complex.
     * @return A constant reference to the map from simplex nodes to integer values.
     */
    const node_stack& get_stack() const;

     /**
     * @brief Update the node_stack by remplacing it.
     * @param new_stack New node_stack to replace the actual one.
     */
    void update_stack(node_stack new_F);


private:
    const SimplexTree& simplex_tree;  ///< Reference to the simplicial complex given in input
    node_stack F; // node_stack values associated to the simplicial complex

    /**
     * @brief Internal function to find an eligible simplex in a list, based on a condition.
     * @param T Boolean map marking simplices.
     * @param simplex_list List of simplices.
     * @param order "increasing" or "decreasing" order condition.
     * @param s_ptr Reference simplex.
     * @return Pointer to the selected simplex or nullptr.
     */
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T, const node_list& simplex_list, std::string order, const node_ptr& s_ptr);

    
    /**
     * @brief Internal function to find an eligible simplex based on filtration values.
     * @param T Boolean map marking simplices.
     * @param simplex_list List of simplices.
     * @param s_ptr Reference simplex.
     * @param F Filtration values.
     * @return Pointer to the selected simplex or nullptr.
     */
    node_ptr find_out(const tsl::robin_map<node_ptr, bool>& T, const node_list& simplex_list, node_ptr s_ptr, const node_stack& F);
};

#endif
