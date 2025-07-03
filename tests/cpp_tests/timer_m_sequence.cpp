/*
    Modification of source code st_iterators.hpp, st_filtration.hpp, st.hpp (line 4 - changed import path of simplextree)
*/
#include "../../src/morse_sequence/_core/morse_sequence.h"

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

void write_timing_result(int size, const std::string& operation, const std::string& lang, double time_sec) {
    std::ofstream file("/home/kali/Téléchargements/Morse-Sequence/graphs/timing_results_cpp.csv", std::ios::app);
    static bool header_written = false;
    if (!header_written) {
        file << "size,operation,language,time_seconds\n";
        header_written = true;
    }
    file << size << "," << operation << "," << lang << "," << time_sec << "\n";
}

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */
void timer_m_sequence() {
    std::ofstream("/home/kali/Téléchargements/Morse-Sequence/graphs/timing_results_cpp.csv", std::ios::trunc).close();  // reset CSV

    auto global_start = std::chrono::high_resolution_clock::now();
    printf("timer_m_sequence : start \n\n");

    std::vector<int> list_faces = {10, 20, 50, 60, 75, 100};

    for (int k : list_faces) {
        std::cout << "\n======================= Grid case " << k << " x " << k << " =======================\n\n";

        // === Block: Grid generation ===
        auto block_start = std::chrono::high_resolution_clock::now();

        SimplexList L = MakeFacesVectorized1(k, k);
        SimplexTree st;
        for (simplex_t s : L){
            st.insert(s);
        }
        MorseSequence ms(st);

        auto block_end = std::chrono::high_resolution_clock::now();
        double time_grid = std::chrono::duration<double>(block_end - block_start).count();
        std::cout << "grid creation + MorseSequence C++ : " << time_grid << " seconds\n\n";
        write_timing_result(k, "init", "cpp", time_grid);

        // Decreasing
        auto start_dec = std::chrono::high_resolution_clock::now();
        auto [ms_dec, n_crit_dec] = ms.decreasing(st);
        auto end_dec = std::chrono::high_resolution_clock::now();
        double time_dec = std::chrono::duration<double>(end_dec - start_dec).count();
        std::cout << "decreasing C++ : " << std::fixed << std::setprecision(6) << time_dec << " seconds\n";
        write_timing_result(k, "decreasing", "cpp", time_dec);

        // Increasing
        auto start_crois = std::chrono::high_resolution_clock::now();
        auto [ms_crois, n_crit_crois] = ms.increasing(st);
        auto end_crois = std::chrono::high_resolution_clock::now();
        double time_crois = std::chrono::duration<double>(end_crois - start_crois).count();
        std::cout << "increasing C++ : " << time_crois << " seconds\n\n";
        write_timing_result(k, "increasing", "cpp", time_crois);

        // Max
        block_start = std::chrono::high_resolution_clock::now();

        auto S_max = ms.simplices(std::nullopt);
        std::sort(S_max.begin(), S_max.end(), [&st](const node_ptr& a, const node_ptr& b) {
            return (st.depth(a) < st.depth(b)) || (st.depth(a) == st.depth(b) && st.full_simplex(a) < st.full_simplex(b) );
        });
        tsl::robin_map<node_ptr, int> F_max;
        for (const auto& s : S_max){
            F_max[s] = 0;
        }
        std::sort(S_max.begin(), S_max.end(), [&st, &F_max](const node_ptr& a, const node_ptr& b) {
            return (F_max[a] < F_max[b]) || (F_max[a] == F_max[b] && st.full_simplex(a) < st.full_simplex(b));
        });

        block_end = std::chrono::high_resolution_clock::now();
        double time_max_preprocess = std::chrono::duration<double>(block_end - block_start).count();
        std::cout << "Max preprocessing C++ : " << time_max_preprocess << " seconds\n";
        write_timing_result(k, "max_preprocessing", "cpp", time_max_preprocess);

        auto start_max = std::chrono::high_resolution_clock::now();
        auto [max, n_crit_max] = ms.Max(S_max, F_max);
        auto end_max = std::chrono::high_resolution_clock::now();
        double time_max = std::chrono::duration<double>(end_max - start_max).count();
        std::cout << "max C++ : " << time_max << " seconds\n\n";
        write_timing_result(k, "max", "cpp", time_max);

        // Min
        block_start = std::chrono::high_resolution_clock::now();

        auto S_min = ms.simplices(std::nullopt);
        std::sort(S_min.begin(), S_min.end(), [&st](const node_ptr& a, const node_ptr& b) {
            return (st.depth(a) > st.depth(b)) || (st.depth(a) == st.depth(b) && st.full_simplex(a) > st.full_simplex(b));
        });
        tsl::robin_map<node_ptr, int> F_min;
        for (const auto& s : S_min) F_min[s] = 0;
        std::sort(S_min.begin(), S_min.end(), [&st, &F_min](const node_ptr& a, const node_ptr& b) {
            return (F_min[a] > F_min[b]) || (F_min[a] == F_min[b] && st.full_simplex(a) > st.full_simplex(b));
        });

        block_end = std::chrono::high_resolution_clock::now();
        double time_min_preprocess = std::chrono::duration<double>(block_end - block_start).count();
        std::cout << "Min preprocessing C++ : " << time_min_preprocess << " seconds\n";
        write_timing_result(k, "min_preprocessing", "cpp", time_min_preprocess);

        auto start_min = std::chrono::high_resolution_clock::now();
        auto [min, n_crit_min] = ms.Min(S_min, F_min);
        auto end_min = std::chrono::high_resolution_clock::now();
        double time_min = std::chrono::duration<double>(end_min - start_min).count();
        std::cout << "Min C++ : " << time_min << " seconds\n";
        write_timing_result(k, "min", "cpp", time_min);
    }

    auto global_end = std::chrono::high_resolution_clock::now();
    double time_total = std::chrono::duration<double>(global_end - global_start).count();
    std::cout << "\ntimer_m_sequence : total time = " << time_total << " seconds\n";
    write_timing_result(0, "total", "cpp", time_total);
}




int main() {
    timer_m_sequence();
    return 0;
}

/* Compilation commands (terminal in Morse-Sequence/src/morse_sequence/cpp_tests)
make timer_m_seq

To generate all the test files in one go (terminal in Morse-Sequence/src/morse_sequence/cpp_tests) : 
make 
*/

/*
Github token : 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
Used SSH key instead 
*/
