/*
	Modification du code source st_iterators.hpp, st_filtration.hpp, st.hpp (ligne 4 -import de simplextree modifiée)
	Modification du code source _simplextree.cpp (ligne 1 à 6 + ligne 17 remplacées)
*/

#include "morse_sequence.h"
#include <vector>

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
        if (std::holds_alternative<simplex_t>(item)) {
            auto simplex = std::get<simplex_t>(item);
            std::cout << "Simplex: ";
            for (auto el : simplex) {
                std::cout << el << " ";
            }
            std::cout << std::endl;
        }
        else if (std::holds_alternative<std::pair<simplex_t, simplex_t>>(item)) {
            auto pair = std::get<std::pair<simplex_t, simplex_t>>(item);
            std::cout << "Pair of simplexes: ";
            for (auto el : pair.first) {
                std::cout << el << " ";
            }
            std::cout << "and ";
            for (auto el : pair.second) {
                std::cout << el << " ";
            }
            std::cout << std::endl;
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
        if (std::holds_alternative<simplex_t>(item)) {
            auto simplex = std::get<simplex_t>(item);
            std::cout << "Simplex: ";
            for (auto el : simplex) {
                std::cout << el << " ";
            }
            std::cout << std::endl;
        }
        else if (std::holds_alternative<std::pair<simplex_t, simplex_t>>(item)) {
            auto pair = std::get<std::pair<simplex_t, simplex_t>>(item);
            std::cout << "Pair of simplexes: ";
            for (auto el : pair.first) {
                std::cout << el << " ";
            }
            std::cout << "and ";
            for (auto el : pair.second) {
                std::cout << el << " ";
            }
            std::cout << std::endl;
        }
    }

	printf("\n\n\n");

	printf("Fin du main\n");
	return 0;
}

/* Code de compilation (terminal dans Morse_Frame)
g++ -o test_morse test2.cpp morse_sequence.cpp simplextree-py/simplextree/_simplextree.cpp -std=c++17 -O3 -Wall -lpython3.11
*/

