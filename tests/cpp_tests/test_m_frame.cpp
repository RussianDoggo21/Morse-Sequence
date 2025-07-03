#include "../../src/morse_sequence/_core/morse_sequence.h"

void test_m_frame(){

    printf("test_m_frame : start \n\n");

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


    node_index_map crit_index = ms.generate_critical_index_map(W);
    m_frame reference_map = ms.reference_map(W, crit_index);
    m_frame coreference_map = ms.coreference_map(W, crit_index);

    std::cout << "Reference Map:\n";
    ms.print_m_frame(reference_map, W, crit_index);

    std::cout << "Coreference Map:\n";
    ms.print_m_frame(coreference_map, W, crit_index);
   
}

int main() {
    test_m_frame();
    return 0;
}

/* Compilation commands (terminal in Morse-Sequence/src/morse_sequence/cpp_tests)
make timer_m_seq test_m_seq test_m_frm

To generate all the test files in one go (terminal in Morse-Sequence/src/morse_sequence/cpp_tests) : 
make 
*/