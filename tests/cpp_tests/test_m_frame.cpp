#include "../../src/morse_sequence/_core/morse_sequence.h"
#include <chrono>

void test_m_frame() {

    using clock = std::chrono::high_resolution_clock;
    using duration = std::chrono::duration<double, std::milli>;

    printf("test_m_frame : start \n\n");

    SimplexTree st;

    SimplexList L = {
        {1, 5, 7}, {1, 2, 7},
        {2, 7, 9}, {2, 3, 9},
        {3, 5, 9}, {1, 3, 5},
        {5, 4, 6}, {5, 6, 7},
        {7, 6, 8}, {7, 8, 9},
        {9, 8, 4}, {9, 4, 5},
        {1, 2, 4}, {2, 4, 6},
        {2, 3, 6}, {3, 6, 8},
        {1, 3, 8}, {1, 4, 8}
    };

    for (simplex_t s : L) {
        st.insert(s);
    }

    MorseSequence ms(st);

    auto result = ms.increasing(st);
    m_sequence W = result.first;

    //ms.print_morse_sequence(result, true);

    /*
    std::cout << "\n== Bitmap implementation ==\n";
    auto t_bitmap_total_start = clock::now();

    auto t0 = clock::now();
    node_index_map crit_index = ms.generate_critical_index_map(W);
    auto t1 = clock::now();
    std::cout << "generate_critical_index_map: " << duration(t1 - t0).count() << " ms\n";

    t0 = clock::now();
    m_frame reference_map = ms.reference_map(W, crit_index);
    t1 = clock::now();
    std::cout << "reference_map: " << duration(t1 - t0).count() << " ms\n";

    t0 = clock::now();
    m_frame coreference_map = ms.coreference_map(W, crit_index);
    t1 = clock::now();
    std::cout << "coreference_map: " << duration(t1 - t0).count() << " ms\n";

    std::cout << "\nReference Map (bitmap):\n";
    ms.print_m_frame(reference_map, W, crit_index);

    std::cout << "\nCoreference Map (bitmap):\n";
    ms.print_m_frame(coreference_map, W, crit_index);

    
    auto t_bitmap_total_end = clock::now();
    auto duration_bitmap = t_bitmap_total_end - t_bitmap_total_start;
    std::cout << "Total Bitmap section: " << duration(duration_bitmap).count() << " ms\n";
    */

    /*
    std::cout << "\n== Simplex implementation ==\n";
    auto t_simplex_total_start = clock::now();

    t0 = clock::now();
    m_frame0 ref2 = ms.reference_map0(W);
    t1 = clock::now();
    std::cout << "reference_map0: " << duration(t1 - t0).count() << " ms\n";

    t0 = clock::now();
    m_frame0 coref2 = ms.coreference_map0(W);
    t1 = clock::now();
    std::cout << "coreference_map0: " << duration(t1 - t0).count() << " ms\n";
    
    auto t_simplex_total_end = clock::now();
    auto duration_simplex = t_simplex_total_end - t_simplex_total_start;

    std::cout << "Total Simplex section: " << duration(duration_simplex).count() << " ms\n";
    
    std::cout << "Difference : " << duration(duration_simplex - duration_bitmap).count() << " ms\n";
    */

    m_frame0 ref2 = ms.reference_map0(W);
    std::cout << "\nReference Map (simplex):\n";
    ms.print_m_frame0(ref2, W);

    /*
    node_index_map crit_index = ms.generate_critical_index_map(W);
    m_frame reference_map = ms.reference_map(W, crit_index);
    std::cout << "\nReference Map (bitmap):\n";
    ms.print_m_frame(reference_map, W, crit_index);
    */

}

int main() {
    test_m_frame();
    return 0;
}
