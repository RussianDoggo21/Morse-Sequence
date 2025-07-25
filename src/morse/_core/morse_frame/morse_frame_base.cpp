#include "morse_frame_base.h"

MorseFrameBase::MorseFrameBase(MorseSequence& ms, const m_sequence& W): UnionFindMF(), ms(ms), simplex_tree(ms.get_simplex_tree()) {

    // Extract critical simplices from W
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr crit = std::get<node_ptr>(item);
            critics.push_back(crit);
        }
    }

    // Build critToIndex
    size_t index = 0;
    for (const node_ptr& crit : critics) {
        critToIndex[crit] = index++;
    }

    // Build indexToCrit
    for (const auto& [crit, idx] : critToIndex) {
        indexToCrit[idx] = crit;
    }

    // Set dimension
    dim_crit = critics.size();
}


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

void MorseFrameBase::print_m_frame(const m_sequence& W) const {

    tsl::robin_map<node_ptr, bitmap> bitarray = get_bitarray();
    for (const auto& item : W) {
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr face_ptr = std::get<node_ptr>(item);
            if (face_ptr) {
                std::cout << "Key (Critical simplex): ";
                simplex_tree.print_simplex(std::cout, face_ptr, false);
                std::cout << " -> Value: ";
                print_bitmap(bitarray.at(face_ptr), W);
                std::cout << "\n";
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }

        else if (std::holds_alternative<node_pair>(item)) {
            node_pair pair = std::get<node_pair>(item);

            if (pair.first && pair.second) {
                std::cout << "Pair of simplices:\n";

                std::cout << "Key (Lower pair): ";
                simplex_tree.print_simplex(std::cout, pair.first, false);
                std::cout << " -> Value: ";
                print_bitmap(bitarray.at(pair.first), W);
                std::cout << "\n";

                std::cout << "Key (Upper pair): ";
                simplex_tree.print_simplex(std::cout, pair.second, false);
                std::cout << " -> Value: ";
                print_bitmap(bitarray.at(pair.second), W);
                std::cout << "\n";
            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }

        std::cout << "\n";
    }
}


