/*
	Modification du code source st_iterators.hpp, st_filtration.hpp, st.hpp (ligne 4 -import de simplextree modifiée)
	Modification du code source _simplextree.cpp (ligne 1 à 6 + ligne 17 remplacées)
*/

#include "morse_sequence.h"

using SimplexList = std::vector<simplex_t>;  // Vecteur de simplexes

int main() {
	printf("Début du main\n");
	SimplexTree st;  // Création d'un complexe simplicial
	SimplexList L = {
        {1, 5, 7}, {1, 2, 7},  // Haut gauche
        {2, 7, 9}, {2, 3, 9},  // Haut milieu
        {3, 5, 9}, {1, 3, 5},  // Haut droit
        {5, 4, 6}, {5, 6, 7},  // Milieu gauche
        {7, 6, 8}, {7, 8, 9},  // Milieu centre
        {9, 8, 4}, {9, 4, 5},  // Milieu droit
        {1, 2, 4}, {2, 4, 6},  // Bas gauche
        {2, 3, 6}, {3, 6, 8},  // Bas milieu
        {1, 3, 8}, {1, 4, 8}   // Bas droit
    	};
	
	for (simplex_t s : L)
	{
		st.insert(s);
	}

    /*
    // Trouver les nœuds dans l'arbre
    node_ptr v0 = st.find(simplex_t {1});
    auto v1 = st.find(simplex_t {3,1});
    auto v2 = st.find(simplex_t {3,1,8});

    st.print_simplex(std::cout, v0, true);
    st.print_simplex(std::cout, v1, true);
    st.print_simplex(std::cout, v2, true);
    */

	MorseSequence ms(st); 
	
	printf("\n\n\n");
	
	printf("Séquence de Morse décroissante\n\n");
	auto result = ms.morse_seq_dec(st);
    
	// 3. Extraire les résultats de la fonction
	auto& morse_sequence = result.first;  // Vecteur des simplexes et paires
	int n_crit = result.second;  // Le critère n_crit

	// Afficher les résultats
	std::cout << "Number of critical points: " << n_crit << std::endl;
	
	for (const auto& item : morse_sequence) {
        // Vérifier le type de l'élément
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr face_ptr = std::get<node_ptr>(item);
            // Vérifie ici si face_ptr n'est pas un pointeur nul avant de l'utiliser
            if (face_ptr) {
                std::cout << "Critical simplex : ";
                st.print_simplex(std::cout, face_ptr, true);
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }
        else if (std::holds_alternative<std::pair<node_ptr, node_ptr>>(item)) {
            std::pair<node_ptr, node_ptr> pair = std::get<std::pair<node_ptr, node_ptr>>(item);
            
            if (pair.first && pair.second) {
                std::cout << "Pair of simplexes: ";
                st.print_simplex(std::cout, pair.first, false);
                std::cout <<  " and ";
                st.print_simplex(std::cout, pair.second, true);
                std::cout << std::endl;
            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }
        
    }
    
    printf("\n\n\n");
 	
	printf("Séquence de Morse croissante\n\n");
	auto result2 = ms.morse_seq_crois(st);

    
	// 3. Extraire les résultats de la fonction
	auto& morse_sequence2 = result2.first;  // Vecteur des simplexes et paires
	int n_crit2 = result2.second;  // Le critère n_crit

	// Afficher les résultats
	std::cout << "Number of critical points: " << n_crit2 << std::endl;
	
	
	for (const auto& item : morse_sequence2) {
        // Vérifier le type de l'élément
        if (std::holds_alternative<node_ptr>(item)) {
            node_ptr face_ptr = std::get<node_ptr>(item);
            // Vérifie ici si face_ptr n'est pas un pointeur nul avant de l'utiliser
            if (face_ptr) {
                std::cout << "Critical simplex: ";
                st.print_simplex(std::cout, face_ptr, true);
            } else {
                std::cout << "Null pointer encountered!" << std::endl;
            }
        }
        else if (std::holds_alternative<std::pair<node_ptr, node_ptr>>(item)) {
            std::pair<node_ptr, node_ptr> pair = std::get<std::pair<node_ptr, node_ptr>>(item);
            
            if (pair.first && pair.second) {
                std::cout << "Pair of simplexes: ";
                st.print_simplex(std::cout, pair.first, false);
                std::cout << " and ";
                st.print_simplex(std::cout, pair.second, true);
                std::cout << std::endl;
            } else {
                std::cout << "Null pointer in pair!" << std::endl;
            }
        }
        
    }

	printf("\n\n\n");
    
	printf("Fin du main\n");
	return 0;
}

/*
#include <iostream>
#include "morse_sequence.h"

int main() {
    SimplexTree st;

    // Insérer des simplexes
    st.insert(simplex_t {0});
    st.insert(simplex_t {0,1});
    st.insert(simplex_t {0,1,2});
    st.insert(simplex_t {3,1,2});

    MorseSequence ms(st);

    // Trouver les nœuds dans l'arbre
    auto v0 = st.find(simplex_t {0});
    auto v1 = st.find(simplex_t {0,1});
    auto v2 = st.find(simplex_t {0,1,2});

    // Affichage du boundary
    std::cout << "Affichage du bord de {0, 1, 2}" << std::endl;
    vector<node_ptr> boundary = ms.boundary(v2);
    for (node_ptr face_ptr : boundary){
        st.print_simplex(std::cout, face_ptr, true);
    }

    //Affichage du coboundary
    std::cout << "Affichage du cobord de {0}" << std::endl;
    vector<node_ptr> coboundary = ms.coboundary(v0);
    for (node_ptr coface_ptr : coboundary) {
        st.print_simplex(std::cout, coface_ptr, true);
    }

    return 0;
}
*/
/* Code de compilation (terminal dans Morse_Frame)
g++ -o test_morse main.cpp morse_sequence.cpp simplextree-py/simplextree/_simplextree.cpp -std=c++20 -O3 -Wall -lpython3.11
g++ -o test_restruct main.cpp restructuration.cpp -std=c++20 -O3 -Wall -lpython3.11
Github token : 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
*/

/*
    // Tester la fonction depth
    std::cout << "Depth de {0} : " << st.depth(v0) << std::endl;  // 1
    std::cout << "Depth de {0,1} : " << st.depth(v1) << std::endl; // 2
    std::cout << "Depth de {0,1,2} : " << st.depth(v2) << std::endl; // 3

    // Affichage du simplexe
    st.print_simplex(std::cout, v2, true); // Devrait afficher "{ 0 1 2 }"


    // Degré du simplexe
    // A quel simplexe est lié un id ?
    for (idx_t i = 0; i <= 3; i++){
        auto cn = st.find_by_id(st.root->children, i);
        std::cout << "Print du simplex lié à l'index " << i << std::endl;
        st.print_simplex(std::cout, cn, true);
    }
    

    std::cout << "Degré de 0: " << st.degree(0) << std::endl;
    std::cout << "Degré de 1: " << st.degree(1) << std::endl;
    std::cout << "Degré de 2: " << st.degree(2) << std::endl;
    std::cout << "Degré de 3: " << st.degree(3) << std::endl;
    std::cout << "Degré de ?: " << st.degree(4) << std::endl;
*/