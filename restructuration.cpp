// Idée de morse_sequence : on donne un complexe en entrée
// En sortie : une séquence de Morse

#include "morse_sequence.h"
#include "simplextree-py/include/utility/combinations.h"
//#include "simplextree-py/include/simplextree/st.hpp"
#include <iostream>
#include <stdexcept> // Pour std::invalid_argument
#include <sstream> // Pour std::ostringstream
#include <algorithm> // Pour std::find
#include <optional>
#include <unordered_map>
#include <functional> // Pour std::hash
#include <variant>
#include <unordered_set>

using std::vector;
using std::set;

// Définition du hash pour simplex_t (std::vector<size_t>) (voir fonction principale std::unordered_map)
namespace std 
{
	template <>
	struct hash<std::vector<size_t>> 
	{
		size_t operator()(const std::vector<size_t>& v) const 
		{
			size_t hash_value = v.size();
			for (size_t x : v) 
			{
				hash_value ^= std::hash<size_t>{}(x) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
			}
			return hash_value;
		}
	};
}


// Constructeur de la classe MorseSequence
// Prend en entrée un SimplexTree
// O(1)
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) {
	std::cout << "MorseSequence créée" << std::endl;
}


// MorseSequence::dim <=> SimplexTree::depth(node_ptr cn) : c.f. utilisation dans main.cpp
// MorseSequence::print_simplex <=> SimplexTree::print_simplex : cf main.cpp


// Renvoie les pointeurs des simplexes composant le bord d'un simplexe sigma
vector<node_ptr> MorseSequence::boundary(node_ptr cn) {   
    vector<node_ptr> boundary; // Vecteur à renvoyer
    faces<> bord(&simplex_tree, cn); // Itérateur de faces du simplexe lié au noeud cn
    // Itération sur les faces
    for (auto& face : bord) {
        node_ptr face_ptr = std::get<0>(face); // On récupère le pointeur de chaque face
        if (face_ptr != cn && simplex_tree.depth(face_ptr) == simplex_tree.depth(cn) - 1){
            //simplex_t s = simplex_tree.full_simplex(face_ptr); 
            boundary.push_back(face_ptr); // On l'ajoute au bord
        }    
    }
    return boundary;
}

// Renvoie les pointeurs des simplexes composant le cobord d'un simplexe sigma
vector<node_ptr> MorseSequence::coboundary(node_ptr cn) {   
    vector<node_ptr> coboundary; // Vecteur à renvoyer
    cofaces<> cobord(&simplex_tree, cn); // Itérateur de cofaces du simplexe lié au noeud cn
    // Itération sur les cofaces
    for (auto& coface : cobord) {
        node_ptr coface_ptr = std::get<0>(coface); // On récupère le pointeur de chaque coface
        if (coface_ptr != cn && simplex_tree.depth(coface_ptr) == simplex_tree.depth(cn) + 1){ // On exclut le simplexe dont on cherche le cobord et on trie avec les bonnes dimensions
            //simplex_t s = simplex_tree.full_simplex(coface_ptr); 
            coboundary.push_back(coface_ptr); // On l'ajoute au cobord
        }    
    }
    return coboundary;
}


// Permet de trouver le simplexe v qui est dans B mais pas dans S
// En notant n = S.size() et m = B.size() : O(nm)
std::optional<node_ptr> MorseSequence::find_out(vector<node_ptr> B, set<node_ptr> S) {
	
    vector<node_ptr> possibilities;

	for (node_ptr v : B) // O(m)
	{
		// Vérifie si v n'est pas dans S
		//if (std::find(S.begin(), S.end(), v) == S.end()) // O(n)
        if (!S.contains(v))
		{
			possibilities.push_back(v); // O(1)
		}
	}

	if (possibilities.empty()) 
	{
		return std::nullopt;  // Retourne un "None" en C++
	}

	return possibilities.at(0); // Retourne le premier élément trouvé
}


// Pour trier une liste de simplexe par dimensions décroissantes
// En notant n = L.size() : O(nlog(n))
vector<node_ptr> MorseSequence::tri_dim_decroissant(vector<node_ptr> L) 
{
	std::sort(L.begin(), L.end(), [this](node_ptr a_ptr, node_ptr b_ptr) // O(nlog(n)) : complexité de std::sort
	{
		return simplex_tree.depth(a_ptr) > simplex_tree.depth(b_ptr); // Tri décroissant par dimension
	});

	return L;
}

// Pour trier une liste de simplexe par dimensions croissantes
// En notant n = L.size() : O(nlog(n))
vector<node_ptr> MorseSequence::tri_dim_croissant(vector<node_ptr> L) 
{
	std::sort(L.begin(), L.end(), [this](node_ptr a_ptr, node_ptr b_ptr) // O(nlog(n)) : complexité de std::sort
	{
		return simplex_tree.depth(a_ptr) < simplex_tree.depth(b_ptr); // Tri croissant par dimension
	});

	return L;
}


// Renvoie les p-simplexes du complexe simplicial si p est précisé
// Renvoie tous les simplexes si p n'est pas précisé
// En notant n = simplextree.size() : O(n)
vector<node_ptr> MorseSequence::simplices(std::optional<int> p) const {
	vector<node_ptr> F;

    // La fonction traverse requiert que lambda_f possède trois arguments de type node_ptr, idx_t et const simplex_t&
	auto lambda_f = [&F](node_ptr cn, idx_t, const simplex_t&) -> bool{
        F.push_back(cn);
        return true;
    };

	// Cas où `p` est None → Traverser tout le complexe
	if (!p) {
		auto tr = st::k_skeleton<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), simplex_tree.dimension()); // O(n) : complexité de st::k_skeleton
		traverse(tr, lambda_f); // O(n) =  O(1) * n : complexité de lambda_f * nombres de simplexes 
	}
	
	// Cas où `p` est défini → Traverser uniquement les `p`-simplexes
	else {
		auto tr = st::k_simplices<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), *p); // O(m) avec m < n
		traverse(tr, lambda_f);  // O(m) =  O(1) * m : complexité de lambda_f * nombres de simplexes 
	}

	return F;
}



// Création d'une séquence de Morse décroissante
// En notant n = simplextree.size(), m = la taille d'un simplexe : 
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::morse_seq_dec(const SimplexTree& st){
    // Récupération des simplexes triés par dimension décroissante
    vector<node_ptr> K = this->simplices(); // O(n)



    // Besoin de trier en python aussi
    this->tri_dim_decroissant(K); // O(nlog(n)) 
    

    // Déclaration et initialisation des variables
    size_t n = K.size();
    size_t i = 0;
    int n_crit = 0;

    std::unordered_map<node_ptr, int> ro;
    vector<node_ptr> L;
    set<node_ptr> S;
    vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> W;

    // Calcul du nombre de cobords pour chaque simplex
    // O(nmlog(m))
    for (node_ptr s : K) { // O(n)
        ro[s] = (this->coboundary(s)).size(); // O(mlog(m))
    }

    // Initialisation de L avec les simplexes ayant un seul cobord
    // O(n)
    for (node_ptr s : K){ // O(n)
        if (ro[s] == 1){ // O(1)
            L.push_back(s); // O(1)
        }
    }

    // Boucle principale
    while (i < n) // O(n)
    {
        // Traitement des éléments de L
        while (!L.empty()) // O(n)
        {
            node_ptr tau = L.back(); // O(?)

            L.pop_back(); // O(?)

            std::optional<node_ptr> v_opt = this->find_out(this->coboundary(tau), S); // O(m²log(m))
            if (!v_opt) // Vérifie si v_opt est nullopt (aucun élément trouvé)
            {
                //std::cout << "Aucun élément trouvé pour tau, on continue.\n";
                continue;
            }
            
            node_ptr v = *v_opt; // Déballer la valeur

            W.push_back(std::make_pair(v, tau));

            S.insert(v);
            S.insert(tau);            

            // Construction de la liste des bords
            vector<node_ptr> bords;
            for (node_ptr elmt : this->boundary(v))
            {
                bords.push_back(elmt);
            }
            for (node_ptr elmt : this->boundary(tau))
            {
                bords.push_back(elmt);
            }

            // Mise à jour de ro et ajout des éléments à L si nécessaire
            for (node_ptr mu : bords) 
            {
                ro[mu] -= 1;
                if (ro[mu] == 1) 
                {
                    L.push_back(mu);
                }
            }
        }

        // Ignorer les éléments déjà dans S
        while (i < n && std::find(S.begin(), S.end(), K[i]) != S.end()) 
        {
            i++;
        }

        if (i == n) 
        {
            return {W, n_crit};
        }

        // Ajout d'un simplex critique
        node_ptr critical = K[i];
        W.push_back(critical);
        n_crit++;

        S.insert(critical);
        

        // Mise à jour de ro et ajout des éléments à L si nécessaire
        for (node_ptr simplex : this->boundary(critical)) 
        {
            ro[simplex] -= 1;
            if (ro[simplex] == 1) 
            {
                L.push_back(simplex);
            }
        }

        // Affichage de fin
        if (i == n - 1) 
        {
            std::cout << "len(S) == n\n";
            std::cout << "=================== End of the program ===================\n";
        }
    }
    return {W, n_crit};
}


// Création d'une séquence de Morse croissante
// En notant n = simplextree.size() : 
std::pair<std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>, int> MorseSequence::morse_seq_crois(const SimplexTree& st)
{
    // Trier les simplexes de st par dimension croissante
    vector<node_ptr> K = this->simplices(); // O(n)

    size_t n = K.size();
    size_t i = 0;
    int n_crit = 0;
    
    set<node_ptr> S;
    vector<node_ptr> L;
    vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>> W;
    std::unordered_map<node_ptr, int> ro;
    
    for (const node_ptr& s : K) {
        ro[s] = 0;
    }
    
    while (i < n - 1) {
        node_ptr sigma = K[i];
        S.insert(sigma);
        W.push_back(sigma);
        n_crit++;
        
        for (const node_ptr& tau : this->coboundary(sigma)) {
            ro[tau] += 1;
            L.push_back(tau);
        }
        
        //Besoin de tri en python aussi
        this->tri_dim_decroissant(L);
        
        while (!L.empty()) {

            node_ptr tau = L.back();
            L.pop_back();
            
            if (ro[tau] == static_cast<int>(st.depth(tau) - 1)) { 
                
                std::optional<node_ptr> v_opt = this->find_out(this->boundary(tau), S);

                if (!v_opt.has_value()) continue;
                
                node_ptr v = *v_opt;
                
                W.push_back(std::make_pair(v, tau));
                S.insert(v);
                S.insert(tau);
                
                
                vector<node_ptr> cobords;
                for (const node_ptr& elmt : this->coboundary(v)) {
                    if (std::find(cobords.begin(), cobords.end(), elmt) == cobords.end()) {
                        cobords.push_back(elmt);
                    }
                }
                for (const node_ptr& elmt : this->coboundary(tau)) {
                    if (std::find(cobords.begin(), cobords.end(), elmt) == cobords.end()) {
                        cobords.push_back(elmt);
                    }
                }
                
                for (const node_ptr& mu : cobords) {
                    ro[mu] += 1;
                    if (std::find(L.begin(), L.end(), mu) == L.end()) {
                        L.push_back(mu);
                    }
                }
                
                this->tri_dim_decroissant(L);
                // Besoin de tri en python
                // Tri croissant ou décroissant ??
            }
        }
        
        while (i < n - 1 && S.contains(K[i])) {
            i++;
        }
    }
    return {W, n_crit};
}

