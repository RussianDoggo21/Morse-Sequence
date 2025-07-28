#include "morse_frame/morse_frame_base.h"

/**
 * @brief Constructor initializes critical simplex mapping and bitmap size.
 * @param ms Reference to the MorseSequence.
 * @param W Sequence of critical simplices and pairs.
 */
MorseFrameBase::MorseFrameBase(MorseSequence& ms, const m_sequence& W)
    : UnionFindMF(), ms(ms), simplex_tree(ms.get_simplex_tree()) {

    // Extract critical simplices
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit = std::get<node_ptr>(item);
            critics.push_back(crit);
        }
    }

    // Assign each critical simplex a unique index
    size_t index = 0;
    for (const node_ptr& crit : critics) {
        critToIndex[crit] = index++;
    }

    // Reverse mapping: index to critical simplex
    for (const auto& [crit, idx] : critToIndex) {
        indexToCrit[idx] = crit;
    }

    // Total number of critical simplices
    dim_crit = critics.size();
}

/**
 * @brief Print the bit vector as a list of critical simplices.
 * @param bm Bitmap to print.
 * @param W Sequence of critical simplices and pairs.
 */
void MorseFrameBase::print_bitmap(const bitmap& bm, const m_sequence& W) const {
    bool empty = true;

    for (size_t i = 0; i < bm.size(); ++i) {
        if (bm[i]) {
            empty = false;
            for (const auto& item : W) {
                if (std::holds_alternative<node_ptr>(item)) {
                    node_ptr c = std::get<node_ptr>(item);
                    if (critToIndex.at(c) == i) {
                        simplex_tree.print_simplex(std::cout, c, false);
                        std::cout << " ";
                        break;
                    }
                }
            }
        }
    }

    if (empty) {
        std::cout << "{}";
    }
}

/**
 * @brief Print all entries in the Morse frame with their associated bitmaps.
 * @param W Sequence of critical simplices and pairs.
 */
void MorseFrameBase::print_m_frame(const m_sequence& W) {
    tsl::robin_map<node_ptr, bitmap> bitarray = get_bitarray();

    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr face_ptr = std::get<node_ptr>(item);
            if (face_ptr) {
                node_ptr root = _find(face_ptr);
                std::cout << "Key (Critical simplex): ";
                simplex_tree.print_simplex(std::cout, face_ptr, false);
                std::cout << " -> Value: ";
                print_bitmap(bitarray.at(root), W);
                std::cout << "\n";
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }

        else if (std::holds_alternative<node_pair>(item)) {
            node_pair pair = std::get<node_pair>(item);
            if (pair.first && pair.second) {
                std::cout << "Pair of simplices:\n";

                // Lower simplex
                node_ptr root1 = _find(pair.first);
                std::cout << "Key (Lower pair): ";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " -> Value: ";
                print_bitmap(bitarray.at(root1), W);
                std::cout << "\n";

                // Upper simplex
                node_ptr root2 = _find(pair.second);
                std::cout << "Key (Upper pair): ";
                simplex_tree.print_simplex(std::cout, pair.second, false);
                std::cout << " -> Value: ";
                print_bitmap(bitarray.at(root2), W);
                std::cout << "\n";

            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }

        std::cout << "\n";
    }
}
