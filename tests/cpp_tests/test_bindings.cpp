#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <pybind11/embed.h> 

#include "../../src/morse_sequence/_core/morse_frame.h"
#include "../../src/morse_sequence/_core/bindings.cpp"

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


// Auxiliary function equivalent to pad_and_pack_S_general in Python 
py::array_t<idx_t> pad_and_pack_S_general(const std::vector<simplex_t>& S_list) {
    size_t max_dim = 0;
    for (const auto& s : S_list)
        if (s.size() > max_dim) max_dim = s.size();

    size_t n = S_list.size();
    size_t cols = max_dim;

    py::array_t<idx_t> arr({n, cols});
    auto buf = arr.mutable_unchecked<2>();

    for (size_t i = 0; i < n; ++i) {
        const simplex_t& s = S_list[i];
        for (size_t j = 0; j < cols; ++j) {
            if (j < s.size()) buf(i, j) = s[j];
            else buf(i, j) = -1;  // padding
        }
    }
    return arr;
}

// Auxiliary function equivalent to pad_and_pack_F_general in Python 
py::array_t<idx_t> pad_and_pack_F_general(const std::unordered_map<simplex_t, int, VectorHash>& F_dict) {
    size_t max_dim = 0;
    for (const auto& [s, val] : F_dict)
        if (s.size() > max_dim) max_dim = s.size();

    size_t n = F_dict.size();
    size_t cols = max_dim + 1;  // last column for the weight

    py::array_t<idx_t> arr({n, cols});
    auto buf = arr.mutable_unchecked<2>();

    size_t i = 0;
    for (const auto& [s, val] : F_dict) {
        for (size_t j = 0; j < max_dim; ++j) {
            if (j < s.size()) buf(i, j) = s[j];
            else buf(i, j) = -1;
        }
        buf(i, cols - 1) = val;
        ++i;
    }
    return arr;
}

int main() {

    py::scoped_interpreter guard{};

    std::vector<int> list_faces = {10, 20, 50, 60, 75};

    for (int k : list_faces) {
        std::cout << "\n======================= Grid case " << k << " x " << k << " =======================\n\n";

        SimplexList L = MakeFacesVectorized1(k, k);

        SimplexTree st;
        for (const auto& s : L) st.insert(s);

        MorseSequence ms(st);

        // --- Préparer S et F pour Max ---

        // Liste simplexes triée par (longueur, lexicographique)
        std::vector<simplex_t> S_max1;
        for (const node_ptr& cn : ms.simplices(std::nullopt)){
            S_max1.push_back(st.full_simplex(cn));
        }
        std::sort(S_max1.begin(), S_max1.end(), [](const simplex_t& a, const simplex_t& b) {
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });

        // Initialiser poids = 0 pour tous
        std::unordered_map<simplex_t, int, VectorHash> F_max1;
        for (const auto& s : S_max1) F_max1[s] = 0;

        // Trier S_max1 par (F_max1[s], taille)
        std::sort(S_max1.begin(), S_max1.end(), [&](const simplex_t& a, const simplex_t& b) {
            if (F_max1[a] != F_max1[b]) return F_max1[a] < F_max1[b];
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });

        // Construire buffers numpy padés
        py::array_t<idx_t> S_buffer_max = pad_and_pack_S_general(S_max1);
        py::array_t<idx_t> F_buffer_max = pad_and_pack_F_general(F_max1);

        // Benchmark _Max_buffered
        auto t0_max = high_resolution_clock::now();
        py::tuple result_max = _Max_buffered(ms, S_buffer_max, F_buffer_max);
        auto t1_max = high_resolution_clock::now();
        std::cout << "_Max_buffered time: " << duration<double>(t1_max - t0_max).count() << " s\n";


        // --- Préparer S et F pour Min ---

        std::vector<simplex_t> S_min1;
        for (const node_ptr& cn : ms.simplices(std::nullopt)){
            S_min1.push_back(st.full_simplex(cn));
        }
        // Tri initial (len, lex) inversé
        std::sort(S_min1.begin(), S_min1.end(), [](const simplex_t& a, const simplex_t& b) {
            if (a.size() != b.size()) return a.size() < b.size();
            return a < b;
        });
        std::reverse(S_min1.begin(), S_min1.end());

        std::unordered_map<simplex_t, int, VectorHash> F_min1;
        for (const auto& s : S_min1) F_min1[s] = 0;

        // Trier S_min1 par (-F_min1[s], -len)
        std::sort(S_min1.begin(), S_min1.end(), [&](const simplex_t& a, const simplex_t& b) {
            if (F_min1[a] != F_min1[b]) return F_min1[a] > F_min1[b];
            if (a.size() != b.size()) return a.size() > b.size();
            return a > b;
        });

        py::array_t<idx_t> S_buffer_min = pad_and_pack_S_general(S_min1);
        py::array_t<idx_t> F_buffer_min = pad_and_pack_F_general(F_min1);

        // Benchmark _Min_buffered
        auto t0_min = high_resolution_clock::now();
        py::tuple result_min = _Min_buffered(ms, S_buffer_min, F_buffer_min);
        auto t1_min = high_resolution_clock::now();
        std::cout << "_Min_buffered time: " << duration<double>(t1_min - t0_min).count() << " s\n";
    }

    return 0;
}

/* Compilation commands (terminal in Morse-Sequence/src/morse_sequence/cpp_tests)
make test_bind

To generate all the test files in one go (terminal in Morse-Sequence/src/morse_sequence/cpp_tests) : 
make 
*/