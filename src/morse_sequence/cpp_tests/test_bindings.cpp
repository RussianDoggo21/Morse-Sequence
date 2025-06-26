#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <pybind11/embed.h> 
#include <chrono>
#include <iostream>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "../_core/morse_sequence.h"
#include "../_core/bindings.cpp"

using simplex_t = SimplexTree::simplex_t;
using node_ptr = SimplexTree::node*;
using m_sequence = std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>;
using node_list = std::vector<node_ptr>;
using m_frame = std::unordered_map<node_ptr, node_list>;

namespace py = pybind11;
using namespace std::chrono;


struct VectorHash {
    std::size_t operator()(const simplex_t& v) const {
        std::size_t seed = v.size();
        for (auto& i : v)
            seed ^= std::hash<idx_t>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

SimplexList MakeFacesVectorized1(int Nr, int Nc) {
    SimplexList out;
    out.reserve(2 * (Nr - 1) * (Nc - 1));  // We know in advance the number of triangles

    auto r = [=](int i, int j) -> unsigned long { 
        return static_cast<unsigned long>(i * Nc + j); 
    };

    for (int i = 0; i < Nr - 1; ++i) {
        for (int j = 0; j < Nc - 1; ++j) {
            // First triangle
            simplex_t tri1 = {r(i, j), r(i, j + 1), r(i + 1, j)};
            // Second triangle
            simplex_t tri2 = {r(i + 1, j), r(i, j + 1), r(i + 1, j + 1)};
            out.push_back(tri1);
            out.push_back(tri2);
        }
    }

    return out;
}


int main() {

    py::scoped_interpreter guard{};

    std::vector<int> list_faces = {10, 20, 50, 60, 75};

    for (int k : list_faces) {
        std::cout << "\n======================= Grid case " << k << " x " << k << " =======================\n\n";

        SimplexList L = MakeFacesVectorized1(k, k);
    

        SimplexTree st;
        /*
        // Remplir le SimplexTree
        std::vector<simplex_t> L = {
            {1, 5, 7}, {1, 2, 7}, {2, 7, 9}, {2, 3, 9}, {3, 5, 9}, {1, 3, 5},
            {5, 4, 6}, {5, 6, 7}, {7, 6, 8}, {7, 8, 9}, {9, 8, 4}, {9, 4, 5},
            {1, 2, 4}, {2, 4, 6}, {2, 3, 6}, {3, 6, 8}, {1, 3, 8}, {1, 4, 8}
        };
        */
        for (const auto& s : L) st.insert(s);

        MorseSequence ms(st);

        // ----------- Partie inspirée du code Python -----------
        std::vector<simplex_t> S_min1;
        for (const node_ptr& sigma : ms.simplices(std::nullopt)) S_min1.push_back(st.full_simplex(sigma));

        // Initialiser F_min1 à zéro
        std::unordered_map<simplex_t, int, VectorHash> F_min1;
        for (const auto& sigma : S_min1) F_min1[sigma] = 0;

        // Trier comme en Python : par valeur (ici toujours 0), puis par -dim, puis lexicographiquement
        std::sort(S_min1.begin(), S_min1.end(), [&](const simplex_t& a, const simplex_t& b) {
            if (F_min1[a] != F_min1[b]) return F_min1[a] > F_min1[b];
            if (a.size() != b.size()) return a.size() > b.size();
            return a > b; // lexicographic
        });

        // Regrouper S et F par dimension
        std::map<int, std::vector<simplex_t>> S_by_dim_min;
        std::map<int, std::vector<std::pair<simplex_t, int>>> F_by_dim_min;

        for (const auto& s : S_min1)
            S_by_dim_min[(int)s.size() - 1].push_back(s);
        for (const auto& [sigma, val] : F_min1)
            F_by_dim_min[(int)sigma.size() - 1].emplace_back(sigma, val);

        // Créer les py::array_t<idx_t>
        std::unordered_map<int, py::array_t<idx_t>> S_arrays_min, F_arrays_min;

        for (const auto& [dim, simplices] : S_by_dim_min) {
            size_t n = simplices.size();
            size_t cols = dim + 1;
            py::array_t<idx_t> arr({n, cols});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n; ++i)
                for (int j = 0; j <= dim; ++j)
                    buf(i, j) = simplices[i][j];
            S_arrays_min[dim] = std::move(arr);
        }

        for (const auto& [dim, pairs] : F_by_dim_min) {
            size_t n = pairs.size();
            size_t cols = dim + 2;  // +1 colonne pour le poids
            py::array_t<idx_t> arr({n, cols});
            auto buf = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n; ++i) {
                const auto& simplex = pairs[i].first;
                int val = pairs[i].second;
                for (int j = 0; j <= dim; ++j)
                    buf(i, j) = simplex[j];
                buf(i, cols - 1) = val;
            }
            F_arrays_min[dim] = std::move(arr);
        }

        // ----------- Benchmark de _Min et _Max (version bindings.cpp) -----------

        auto n0 = high_resolution_clock::now();
        py::tuple result_max = _Max(ms, S_arrays_min, F_arrays_min);
        auto n1 = high_resolution_clock::now();
        std::cout << "_Max time: " << duration<double>(n1 - n0).count() << " s\n";
        
        auto t0 = high_resolution_clock::now();
        py::tuple result_min = _Min(ms, S_arrays_min, F_arrays_min);
        auto t1 = high_resolution_clock::now();
        std::cout << "_Min time: " << duration<double>(t1 - t0).count() << " s\n";


    }

    return 0;
}
