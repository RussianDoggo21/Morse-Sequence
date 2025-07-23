#include "../../src/morse_sequence/_core/morse_frame.h"
#include <chrono>
#include <iomanip>

// Même générateur que timer_m_sequence
SimplexList MakeFacesVectorized1(int Nr, int Nc) {
    SimplexList out;
    out.reserve(2 * (Nr - 1) * (Nc - 1));
    auto r = [=](int i, int j) -> unsigned long {
        return static_cast<unsigned long>(i * Nc + j);
    };

    for (int i = 0; i < Nr - 1; ++i) {
        for (int j = 0; j < Nc - 1; ++j) {
            simplex_t tri1 = {r(i, j), r(i, j + 1), r(i + 1, j)};
            simplex_t tri2 = {r(i + 1, j), r(i, j + 1), r(i + 1, j + 1)};
            out.push_back(tri1);
            out.push_back(tri2);
        }
    }

    return out;
}

void timer_m_frame() {

    using clock = std::chrono::high_resolution_clock;
    using duration = std::chrono::duration<double, std::milli>;

    printf("test_m_frame : start \n\n");

    std::vector<int> grid_sizes = {100, 150, 200, 255};

    for (int k : grid_sizes) {
        std::cout << "\n========== Grid " << k << " x " << k << " ==========\n";

        // Grid + MorseSequence
        auto t_start = clock::now();
        SimplexList L = MakeFacesVectorized1(k, k);
        SimplexTree st;
        for (simplex_t s : L) st.insert(s);
        MorseSequence ms(st);
        auto t_end = clock::now();
        std::cout << "Grid + MorseSequence construction: " << duration(t_end - t_start).count() << " ms\n";

        // Increasing
        t_start = clock::now();
        auto [W, n_crit] = ms.increasing(st);
        t_end = clock::now();
        std::cout << "increasing(): " << duration(t_end - t_start).count() << " ms\n";

        MorseFrame mf(ms);
        
        std::cout << "\n== Simplex implementation ==\n";
        auto t_simplex_total_start = clock::now();

        auto m0 = clock::now();
        m_frame0 ref2 = mf.reference_map0(W);
        auto m1 = clock::now();
        std::cout << "reference_map0: " << duration(m1 - m0).count() << " ms\n";

        m0 = clock::now();
        m_frame0 coref2 = mf.coreference_map0(W);
        m1 = clock::now();
        std::cout << "coreference_map0: " << duration(m1 - m0).count() << " ms\n";

        auto t_simplex_total_end = clock::now();

        auto t_simplex = t_simplex_total_end - t_simplex_total_start;
        std::cout << "Total Simplex section: " << duration(t_simplex).count() << " ms\n";

        std::cout << "\n== Bitmap implementation ==\n";
        auto t_bitmap_total_start = clock::now();

        auto t0 = clock::now();
        node_index_map crit_index = mf.generate_critical_index_map(W);
        auto t1 = clock::now();
        std::cout << "generate_critical_index_map: " << duration(t1 - t0).count() << " ms\n";

        t0 = clock::now();
        m_frame reference_map = mf.reference_map(W, crit_index);
        t1 = clock::now();
        std::cout << "reference_map: " << duration(t1 - t0).count() << " ms\n";

        t0 = clock::now();
        m_frame coreference_map = mf.coreference_map(W, crit_index);
        t1 = clock::now();
        std::cout << "coreference_map: " << duration(t1 - t0).count() << " ms\n";

        auto t_bitmap_total_end = clock::now();

        auto t_bitmap = t_bitmap_total_end - t_bitmap_total_start;
        std::cout << "Total Bitmap section: " << duration(t_bitmap).count() << " ms\n";
        
        std::cout << "Difference: (simplex - bitmap) :" << duration(t_simplex - t_bitmap).count() << " ms\n";
        std::cout << "Ratio (simplex/bitmap): " << static_cast<double>(t_simplex.count()) / t_bitmap.count() << "\n";

                  
    }
}

int main() {
    timer_m_frame();
    return 0;
}
