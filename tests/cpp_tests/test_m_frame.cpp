#include "morse_frame/ref_map.h"
#include "morse_frame/coref_map.h"

void test_m_frame() {

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

    ms.print_morse_sequence(result, true);

    printf("\n ===== Reference Map ===== \n ");
    RefMap ref_map(ms, W);
    ref_map.print_m_frame(W);

    printf("\n ===== Coreference Map ===== \n ");
    CorefMap coref_map(ms, W);
    coref_map.print_m_frame(W);    
}

int main() {
    test_m_frame();
    return 0;
}
