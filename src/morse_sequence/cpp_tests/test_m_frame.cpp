#include "../_core/morse_sequence.h"
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
using m_frame = std::unordered_map<node_ptr, node_list>;

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

    m_frame reference_map = ms.reference_map(W);
    m_frame coreference_map = ms.coreference_map(W);

    std::cout << "Reference Map:\n";
    ms.print_m_frame(reference_map, W);

    std::cout << "Coreference Map:\n";
    ms.print_m_frame(coreference_map, W);
   
}

int main() {
    test_m_frame();
    return 0;
}

/* Compilation commands (terminal in Morse-Sequence/src/morse_sequence)
make timer_m_seq test_m_seq test_m_frm

To generate all the test files in one go (terminal in MOrse-Sequence/src/morse_sequence) : 
make 
*/