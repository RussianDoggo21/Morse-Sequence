#include "morse_frame/ref_map.h"
#include "morse_frame/coref_map.h"

void test_m_frame() {

    printf("test_m_frame : start \n\n");

    SimplexTree st;

    /*
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
    */

    SimplexList L = {
        {1,3,5}, {1,5,6}, {1,3,6},
        {2,3,5}, {2,4,5},
        {1,2,4}, {1,3,4},
        {2,3,8}, {3,4,8},
        {1,2,8}, {1,7,8},
        {1,2,7}, {2,3,7},
        {4,6,7}, {4,7,8},
        {3,6,7}, {4,5,6}
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

    printf("\n ===== Persistence ===== \n ");
    RefMap morse_frame1 = ref_map;
    auto result_ref = morse_frame1.persistence();
    morse_frame1.print_persistence_results(result_ref);

    printf("\n ===== Coreference Map ===== \n ");
    CorefMap coref_map(ms, W);
    coref_map.print_m_frame(W);    

    printf("\n ===== Copersistence ===== \n ");
    CorefMap morse_frame2 = coref_map;
    auto result_coref = morse_frame2.copersistence();
    morse_frame2.print_persistence_results(result_coref);
}

int main() {
    test_m_frame();
    return 0;
}
