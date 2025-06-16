/*
    Modification of source code st_iterators.hpp, st_filtration.hpp, st.hpp (line 4 - changed import path of simplextree)
*/

#include "../morse_sequence.h"
#include <vector>
#include <unordered_map>
#include <algorithm> 
#include <tuple>
#include <cstddef>
#include <functional>
#include <chrono>
#include <iomanip>
using SimplexList = std::vector<simplex_t>;  // Vector of simplices
using m_sequence = std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>;
using node_list = std::vector<node_ptr>;
using morse_frame = std::unordered_map<node_ptr, node_list>;

/*
// Used in order to create types of the form unordered_map<simplex_t, xyz>
// Done by Chat GPT
namespace std {
    template<>
    struct hash<std::vector<std::size_t>> {
        std::size_t operator()(const std::vector<std::size_t>& v) const {
            std::size_t seed = v.size();
            for (auto& i : v) {
                seed ^= std::hash<std::size_t>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
*/

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

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */
void timer_comparison() {
    std::vector<int> list_faces = {10, 20, 50, 60, 75};

    for (int k : list_faces) {
        std::cout << "\n======================= Grid case " << k << " x " << k << " =======================\n\n";

        SimplexList L = MakeFacesVectorized1(k, k);
        SimplexTree st;
        for (simplex_t s : L){
            st.insert(s);
        }
        MorseSequence ms(st);

        // Decreasing
        auto start_dec = std::chrono::high_resolution_clock::now();
        auto [ms_dec, n_crit_dec] = ms.decreasing(st);
        auto end_dec = std::chrono::high_resolution_clock::now();
        double time_dec = std::chrono::duration<double>(end_dec - start_dec).count();
        std::cout << "decreasing C++ : " << std::fixed << std::setprecision(6) << time_dec << " seconds\n";

        // Increasing
        auto start_crois = std::chrono::high_resolution_clock::now();
        auto [ms_crois, n_crit_crois] = ms.increasing(st);
        auto end_crois = std::chrono::high_resolution_clock::now();
        double time_crois = std::chrono::duration<double>(end_crois - start_crois).count();
        std::cout << "increasing C++ : " << time_crois << " seconds\n";

        // Max
        auto S_max = ms.simplices(std::nullopt);
        std::sort(S_max.begin(), S_max.end(), [&st](const node_ptr& a, const node_ptr& b) {
            return (st.depth(a) < st.depth(b)) || (st.depth(a) == st.depth(b) && st.full_simplex(a) < st.full_simplex(b) );
        });
        std::unordered_map<node_ptr, int> F_max;
        for (const auto& s : S_max){
            F_max[s] = 0;
        } 
        std::sort(S_max.begin(), S_max.end(), [&st, &F_max](const node_ptr& a, const node_ptr& b) {
            return (F_max[a] < F_max[b]) || (F_max[a] == F_max[b] && st.full_simplex(a) < st.full_simplex(b));
        });

        auto start_max = std::chrono::high_resolution_clock::now();
        auto [max, n_crit_max] = ms.Max(S_max, F_max);
        auto end_max = std::chrono::high_resolution_clock::now();
        double time_max = std::chrono::duration<double>(end_max - start_max).count();
        std::cout << "max C++ : " << time_max << " seconds\n";

        // Min
        auto S_min = ms.simplices(std::nullopt);
        std::sort(S_min.begin(), S_min.end(), [&st](const node_ptr& a, const node_ptr& b) {
            return (st.depth(a) > st.depth(b)) || (st.depth(a) == st.depth(b) && st.full_simplex(a) > st.full_simplex(b));
        });
        std::unordered_map<node_ptr, int> F_min;
        for (const auto& s : S_min) F_min[s] = 0;
        std::sort(S_min.begin(), S_min.end(), [&st, &F_min](const node_ptr& a, const node_ptr& b) {
            return (F_min[a] > F_min[b]) || (F_min[a] == F_min[b] && st.full_simplex(a) > st.full_simplex(b));
        });

        auto start_min = std::chrono::high_resolution_clock::now();
        auto [min, n_crit_min] = ms.Min(S_min, F_min);
        auto end_min = std::chrono::high_resolution_clock::now();
        double time_min = std::chrono::duration<double>(end_min - start_min).count();
        std::cout << "Min C++ : " << time_min << " seconds\n";
    }
}


void test(){
	SimplexTree st;  // Creation of a simplicial complex
    //SimplexList L = {{0,1,2}};
    //SimplexList L = MakeFacesVectorized1(75, 75);
    
    SimplexList L = {
                        {1, 5, 7}, {1, 2, 7},  // Top left
                        {2, 7, 9}, {2, 3, 9},  // Top middle
                        {3, 5, 9}, {1, 3, 5},  // Top right
                        {5, 4, 6}, {5, 6, 7},  // Middle left
                        {7, 6, 8}, {7, 8, 9},  // Middle center
                        {9, 8, 4}, {9, 4, 5},  // Middle right
                        {1, 2, 4}, {2, 4, 6},  // Bottom left
                        {2, 3, 6}, {3, 6, 8},  // Bottom middle
                        {1, 3, 8}, {1, 4, 8}   // Bottom right
                    };
	for (simplex_t s : L){
		st.insert(s);
	}

	MorseSequence ms(st); 
	
	printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    printf("Decreasing Morse sequence: \n\n");
    auto start_dec = std::chrono::high_resolution_clock::now();
	auto result_dec = ms.decreasing(st);
    auto end_dec = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_dec = end_dec - start_dec;
    std::cout << "Execution time: " << duration_dec.count() << " ms" << std::endl;

    ms.print_morse_sequence(result_dec);
    
    printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    printf("Increasing Morse sequence: \n\n");
    auto start_crois = std::chrono::high_resolution_clock::now();
	auto result_increasing = ms.increasing(st);
    auto end_crois = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_crois = end_crois - start_crois;
    std::cout << "Execution time: " << duration_crois.count() << " ms" << std::endl;
  
    ms.print_morse_sequence(result_increasing);
    
    printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    // S = sorted(ms.simplices(), key=lambda x: (len(x), x))
    node_list S = ms.simplices(std::nullopt);

    std::sort(S.begin(), S.end(), [st](node_ptr a_ptr, node_ptr b_ptr) {
        if (st.depth(a_ptr) != st.depth(b_ptr)) return st.depth(a_ptr) < st.depth(b_ptr);
        return a_ptr < b_ptr; // lexicographical comparison
    });

    // F = dict()
    std::unordered_map<node_ptr, int> F;
    for (node_ptr cn : S) {
       F[cn] = 0;
    }

    // Check if a simplex is included in another => function st.is_face(sigma, tau)

    // S = sorted(S, key=lambda s: (F[s], len(s)))
    std::sort(S.begin(), S.end(), [&F, st](node_ptr a_ptr, node_ptr b_ptr) {
        if (F[a_ptr] != F[b_ptr]) return F[a_ptr] < F[b_ptr];
        return st.depth(a_ptr) < st.depth(b_ptr);
    });

	
    printf("Max: \n\n");
    auto start_max = std::chrono::high_resolution_clock::now();
	auto result_max = ms.Max(S, F);
    auto end_max = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_max = end_max - start_max;
    std::cout << "Execution time: " << duration_max.count() << " ms" << std::endl;
   
    ms.print_morse_sequence(result_max);
    
    printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    // 1. Retrieve simplices and sort by decreasing dimension then lexicographically decreasing
    node_list S2 = ms.simplices(std::nullopt);

    std::sort(S2.begin(), S2.end(), [st](node_ptr a_ptr, node_ptr b_ptr) {
        if (st.depth(a_ptr) != st.depth(b_ptr)) return st.depth(a_ptr) > st.depth(b_ptr); // decreasing order
        return b_ptr < a_ptr; // lexicographically decreasing order
    });

    // 2. Initialize F
    std::unordered_map<node_ptr, int> F2;
    for (node_ptr cn : S2) {
        F2[cn] = 0;
    }

    // 3. Manually set a value of F if needed (example here: F[(0, 1, 2)] = 0)

    // 4. Sort S again according to priority on decreasing F, then decreasing dimension
    std::sort(S2.begin(), S2.end(), [&F2, st](node_ptr a_ptr, node_ptr b_ptr) {
        if (F2[a_ptr] != F2[b_ptr]) return F2[a_ptr] > F2[b_ptr]; // priority on decreasing F
        return st.depth(a_ptr) > st.depth(b_ptr); // then decreasing dimension
    });


	printf("Min\n\n");
    auto start_min = std::chrono::high_resolution_clock::now();
	auto result_min = ms.Min(S2, F2);
    auto end_min = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_min = end_min - start_min;
    std::cout << "Execution time: " << duration_min.count() << " ms" << std::endl;

	ms.print_morse_sequence(result_min);

	printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */
}

void test_reference_map(){
    SimplexTree st;  // Creation of a simplicial complex
    //SimplexList L = {{0,1,2}};

    
    SimplexList L = {
                        {1, 5, 7}, {1, 2, 7},  // Top left
                        {2, 7, 9}, {2, 3, 9},  // Top middle
                        {3, 5, 9}, {1, 3, 5},  // Top right
                        {5, 4, 6}, {5, 6, 7},  // Middle left
                        {7, 6, 8}, {7, 8, 9},  // Middle center
                        {9, 8, 4}, {9, 4, 5},  // Middle right
                        {1, 2, 4}, {2, 4, 6},  // Bottom left
                        {2, 3, 6}, {3, 6, 8},  // Bottom middle
                        {1, 3, 8}, {1, 4, 8}   // Bottom right
                    };
    
	for (simplex_t s : L){
		st.insert(s);
	}

	MorseSequence ms(st); 

    auto result = ms.increasing(st);
    m_sequence W = result.first;

    ms.print_morse_sequence(result, true);


    morse_frame reference_map = ms.reference_map(W);
    morse_frame coreference_map = ms.coreference_map(W);

    std::cout << "Reference Map:\n";
    ms.print_morse_frame(reference_map, W);

    std::cout << "Coreference Map:\n";
    ms.print_morse_frame(coreference_map, W);
   
}

int main() {
	printf("Start of main\n\n");
    test_reference_map();
    //test();
    //timer_comparison();
    printf("End of main");
    return 0;
}

/* Compilation commands (terminal in Morse-Sequence/src/morse_sequence)
g++ -o timer tests/timer.cpp morse_sequence.cpp -std=c++20 -O3 -Wall -lpython3.13

Pour lancer perf :
sudo sysctl kernel.perf_event_paranoid=1
perf record --call-graph=dwarf -F 100 ./timer 100 2000


Github token : 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
Used SSH key instead 
*/
