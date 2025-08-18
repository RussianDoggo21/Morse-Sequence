#include "morse_sequence/morse_sequence.h"

void test_m_sequence(){

    printf("test_m_sequence : start \n\n");

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
	auto result_dec = ms.decreasing();
    auto end_dec = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_dec = end_dec - start_dec;
    std::cout << "Execution time: " << duration_dec.count() << " ms" << std::endl;

    ms.print_morse_sequence(result_dec, true);
    
    printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    printf("Increasing Morse sequence: \n\n");
    auto start_crois = std::chrono::high_resolution_clock::now();
	auto result_increasing = ms.increasing();
    auto end_crois = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_crois = end_crois - start_crois;
    std::cout << "Execution time: " << duration_crois.count() << " ms" << std::endl;
  
    ms.print_morse_sequence(result_increasing, true);
    
    printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    // S = sorted(ms.simplices(), key=lambda x: (len(x), x))
    node_list S = ms.simplices(std::nullopt);

    std::sort(S.begin(), S.end(), [st](node_ptr a_ptr, node_ptr b_ptr) {
        if (st.depth(a_ptr) != st.depth(b_ptr)) return st.depth(a_ptr) < st.depth(b_ptr);
        return a_ptr < b_ptr; // lexicographical comparison
    });

    // F = dict()
    
    tsl::robin_map<node_ptr, int> F = ms.get_stack();
    /*
    for (node_ptr cn : S) {
       F[cn] = 0;
    }
    */

    // Check if a simplex is included in another => function st.is_face(sigma, tau)

    // S = sorted(S, key=lambda s: (F[s], len(s)))
    std::sort(S.begin(), S.end(), [&F, st](node_ptr a_ptr, node_ptr b_ptr) {
        if (F[a_ptr] != F[b_ptr]) return F[a_ptr] < F[b_ptr];
        return st.depth(a_ptr) < st.depth(b_ptr);
    });

	
    printf("Max: \n\n");
    auto start_max = std::chrono::high_resolution_clock::now();
	auto result_max = ms.Max(S);
    auto end_max = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_max = end_max - start_max;
    std::cout << "Execution time: " << duration_max.count() << " ms" << std::endl;
   
    ms.print_morse_sequence(result_max, true);
    
    printf("\n\n\n");

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    // 1. Retrieve simplices and sort by decreasing dimension then lexicographically decreasing
    node_list S2 = ms.simplices(std::nullopt);

    std::sort(S2.begin(), S2.end(), [st](node_ptr a_ptr, node_ptr b_ptr) {
        if (st.depth(a_ptr) != st.depth(b_ptr)) return st.depth(a_ptr) > st.depth(b_ptr); // decreasing order
        return b_ptr < a_ptr; // lexicographically decreasing order
    });

    // 2. Initialize F
    tsl::robin_map<node_ptr, int> F2 = ms.get_stack();
    /*
    for (node_ptr cn : S2) {
        F2[cn] = 0;
    }
    */

    // 3. Manually set a value of F if needed (example here: F[(0, 1, 2)] = 0)

    // 4. Sort S again according to priority on decreasing F, then decreasing dimension
    std::sort(S2.begin(), S2.end(), [&F2, st](node_ptr a_ptr, node_ptr b_ptr) {
        if (F2[a_ptr] != F2[b_ptr]) return F2[a_ptr] > F2[b_ptr]; // priority on decreasing F
        return st.depth(a_ptr) > st.depth(b_ptr); // then decreasing dimension
    });


	printf("Min\n\n");
    auto start_min = std::chrono::high_resolution_clock::now();
	auto result_min = ms.Min(S2);
    auto end_min = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_min = end_min - start_min;
    std::cout << "Execution time: " << duration_min.count() << " ms" << std::endl;

	ms.print_morse_sequence(result_min, true);

	printf("\n\n\n");
}

int main() {
    test_m_sequence();
    return 0;
}
/* Compilation commands (terminal in Morse-Sequence/src/morse_sequence/cpp_tests)
make test_m_seq

To generate all the test files in one go (terminal in Morse-Sequence/src/morse_sequence/cpp_tests) : 
make
*/