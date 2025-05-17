/*
	Modification du code source st_iterators.hpp, st_filtration.hpp, st.hpp (ligne 4 - chemin d'import de simplextree modifiée)
	Modification du code source _simplextree.cpp - lignes 588 et 589 : rajout des fonctions find_node et full_simplex via pybind11
    ............................................ - ligne 541 : rajout de la fonction find_node
    ............................................ - ligne 549 : rajout du type node_ptr via pybind11
*/

#include "morse_sequence.h"
#include <vector>
#include <unordered_map>
#include <algorithm> 
#include <tuple>
#include <cstddef>
#include <functional>
#include <chrono>
#include <iomanip>
using SimplexList = std::vector<simplex_t>;  // Vecteur de simplexes

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

std::vector<simplex_t> MakeFacesVectorized1(int Nr, int Nc) {
    std::vector<simplex_t> out;
    out.reserve(2 * (Nr - 1) * (Nc - 1));  // On sait à l'avance le nombre de triangles

    auto r = [=](int i, int j) -> unsigned long { 
        return static_cast<unsigned long>(i * Nc + j); 
    };

    for (int i = 0; i < Nr - 1; ++i) {
        for (int j = 0; j < Nc - 1; ++j) {
            // Premier triangle
            simplex_t tri1 = {r(i, j), r(i, j + 1), r(i + 1, j)};
            // Deuxième triangle
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
        std::cout << "\n======================= Cas grille " << k << " x " << k << " =======================\n\n";

        SimplexList L = MakeFacesVectorized1(k, k);
        SimplexTree st;
        for (simplex_t s : L){
            st.insert(s);
        }
        MorseSequence ms(st);

        // Décroissante
        auto start_dec = std::chrono::high_resolution_clock::now();
        auto [ms_dec, n_crit_dec] = ms.morse_seq_decrois(st);
        auto end_dec = std::chrono::high_resolution_clock::now();
        double time_dec = std::chrono::duration<double>(end_dec - start_dec).count();
        std::cout << "décroissante C++ : " << std::fixed << std::setprecision(6) << time_dec << " secondes\n";

        // Croissante
        auto start_crois = std::chrono::high_resolution_clock::now();
        auto [ms_crois, n_crit_crois] = ms.morse_seq_crois(st);
        auto end_crois = std::chrono::high_resolution_clock::now();
        double time_crois = std::chrono::duration<double>(end_crois - start_crois).count();
        std::cout << "croissante C++ : " << time_crois << " secondes\n";

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
        std::cout << "max C++ : " << time_max << " secondes\n";

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
        std::cout << "Min C++ : " << time_min << " secondes\n";
    }
}

int main() {
	printf("Début du main\n");
    //timer_comparison();
    
	SimplexTree st;  // Création d'un complexe simplicial
    //SimplexList L = {{0,1,2}};
    //SimplexList L = MakeFacesVectorized1(10, 10);
    
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
    
	for (simplex_t s : L){
		st.insert(s);
	}

	MorseSequence ms(st); 
	
	printf("\n\n\n");
    
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */
/*
    printf("Séquence de Morse croissante : \n\n");
    auto start_crois = std::chrono::high_resolution_clock::now();
	auto result = ms.morse_seq_crois(st);
    auto end_crois = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_crois = end_crois - start_crois;
    std::cout << "Temps d'exécution : " << duration_crois.count() << " ms" << std::endl;
 
    
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
*/
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */
/*
    printf("Séquence de Morse décroissante : \n\n");
    auto start_dec = std::chrono::high_resolution_clock::now();
	auto result_dec = ms.morse_seq_decrois(st);
    auto end_dec = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_dec = end_dec - start_dec;
    std::cout << "Temps d'exécution : " << duration_dec.count() << " ms" << std::endl;

    
	// 3. Extraire les résultats de la fonction
	auto& morse_sequence_dec = result_dec.first;  // Vecteur des simplexes et paires
	int n_crit_dec = result_dec.second;  // Le critère n_crit

	// Afficher les résultats
	std::cout << "Number of critical points: " << n_crit_dec << std::endl;
	
	for (const auto& item : morse_sequence_dec) {
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
*/
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    // S = sorted(ms.simplices(), key=lambda x: (len(x), x))
    std::vector<node_ptr> S = ms.simplices(std::nullopt);

    std::sort(S.begin(), S.end(), [st](node_ptr a_ptr, node_ptr b_ptr) {
        if (st.depth(a_ptr) != st.depth(b_ptr)) return st.depth(a_ptr) < st.depth(b_ptr);
        return a_ptr < b_ptr; // lexicographical comparison
    });

    // F = dict()
    std::unordered_map<node_ptr, int> F;
    for (node_ptr cn : S) {
       F[cn] = 0;
    }

    // Vérifier si un simplexe est inclus dans un autre => fonction st.is_face(sigma, tau) 

    // S = sorted(S, key=lambda s: (F[s], len(s)))
    std::sort(S.begin(), S.end(), [&F, st](node_ptr a_ptr, node_ptr b_ptr) {
        if (F[a_ptr] != F[b_ptr]) return F[a_ptr] < F[b_ptr];
        return st.depth(a_ptr) < st.depth(b_ptr);
    });

	
    printf("Max : \n\n");
    auto start_max = std::chrono::high_resolution_clock::now();
	auto result_max = ms.Max(S, F);
    auto end_max = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_max = end_max - start_max;
    std::cout << "Temps d'exécution : " << duration_max.count() << " ms" << std::endl;
   
    
	// 3. Extraire les résultats de la fonction
	auto& morse_sequence_max = result_max.first;  // Vecteur des simplexes et paires
	int n_crit_max = result_max.second;  // Le critère n_crit

	// Afficher les résultats
	std::cout << "Number of critical points: " << n_crit_max << std::endl;
	
	for (const auto& item : morse_sequence_max) {
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

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

    // 1. Récupérer les simplexes et les trier par ordre décroissant de dimension puis lexicographiquement décroissant
    std::vector<node_ptr> S2 = ms.simplices(std::nullopt);

    std::sort(S2.begin(), S2.end(), [st](node_ptr a_ptr, node_ptr b_ptr) {
        if (st.depth(a_ptr) != st.depth(b_ptr)) return st.depth(a_ptr) > st.depth(b_ptr); // ordre décroissant
        return b_ptr < a_ptr; // ordre lexicographique décroissant
    });

    // 2. Initialiser F
    std::unordered_map<node_ptr, int> F2;
    for (node_ptr cn : S2) {
        F2[cn] = 0;
    }

    // 3. Définir manuellement une valeur de F si besoin (exemple ici : F[(0, 1, 2)] = 0)

    // 4. Trier de nouveau S selon la priorité à F décroissant, puis à dimension décroissante
    std::sort(S2.begin(), S2.end(), [&F2, st](node_ptr a_ptr, node_ptr b_ptr) {
        if (F2[a_ptr] != F2[b_ptr]) return F2[a_ptr] > F2[b_ptr]; // priorité à F décroissant
        return st.depth(a_ptr) > st.depth(b_ptr); // puis dimension décroissante
    });


	printf("Min\n\n");
    auto start_min = std::chrono::high_resolution_clock::now();
	auto result2 = ms.Min(S2, F2);
    auto end_min = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_min = end_min - start_min;
    std::cout << "Temps d'exécution : " << duration_min.count() << " ms" << std::endl;

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

/* ----------------------------------------------------------------------------------------------------------------------------------------------------------- */

	printf("Fin du main\n");
	return 0;
}

/* Code de compilation (terminal dans Morse_Frame)
g++ -o test_morse main.cpp old_morse_sequence.cpp simplextree-py/simplextree/_simplextree.cpp -std=c++20 -O3 -Wall -lpython3.13
g++ -o test_restruct main.cpp restructuration.cpp -std=c++20 -O3 -Wall -lpython3.13
g++ -o f_sequence main.cpp morse_sequence.cpp -std=c++20 -O3 -Wall -lpython3.13

Github token : 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
Utilisé la clé SSH à la place 
*/




