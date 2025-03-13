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
// O(1)
MorseSequence::MorseSequence(const SimplexTree& st) : simplex_tree(st) 
{
	std::cout << "MorseSequence créée" << std::endl;
}

// Renvoie la dimension d'un simplexe sigma
// O(1)
int MorseSequence::dim(const simplex_t& sigma)
{
	return sigma.size() - 1;
}

// Fonction pour print un simplexe 
// En notant n = sigma.size() : O(n)
void MorseSequence::print_simplex(const simplex_t& sigma) const 
{
	std::cout << "{";
	for (size_t i = 0; i < sigma.size(); ++i) 
	{
		std::cout << sigma[i];
		if (i < sigma.size() - 1) 
		{
	    		std::cout << ", ";
		}
	}
	std::cout << "}" << std::endl; // Ajout du saut de ligne
}

// Renvoie le bord d'un simplexe sigma
// En notant n = sigma.size() : O(n²)
vector<simplex_t> MorseSequence::boundary(const simplex_t& sigma) 
{
    vector<simplex_t> result;
    simplex_t sigma_copy = sigma;

    for_each_combination // O(n) : boucle sur n éléments
    (
        sigma_copy.begin(), sigma_copy.begin() + (sigma_copy.size() - 1), sigma_copy.end(),
        [&result](auto first, auto mid) 
        {
            result.emplace_back(first, mid);  // O(n) : copie d'un simplex de taille n
            return false; // Continuer l'itération
        }
    );

    return result;
}

// Renvoie les cofaces d'un simplexe sigma
// En notant n = sigma.size() : O(n.log(n))
vector<simplex_t> MorseSequence::cofaces(const simplex_t& sigma) 
{
	vector<simplex_t> F;

	// Vérifier si le simplexe est vide
	if (sigma.empty()) 
	{
		throw std::invalid_argument("Erreur : sigma ne peut pas être une liste vide.");
	}

	// Vérifier si sigma est dans le complexe simplicial
	if (!simplex_tree.find(sigma)) 
	{
		throw std::invalid_argument("Erreur : sigma n'appartient pas au complexe simplicial.");
	}

	// Lambda pour stocker les cofaces, en excluant sigma
	auto lambda_f = [this, &F, &sigma](node_ptr, idx_t, const simplex_t& s) -> bool
	{
		simplex_t sorted_s = s;
		simplex_t sorted_sigma = sigma;

		std::sort(sorted_s.begin(), sorted_s.end()); // O(n.log(n)) : tri d'un vecteur
		std::sort(sorted_sigma.begin(), sorted_sigma.end()); // O(n.log(n)) : tri d'un vecteur
		
		// On n'ajoute que les cofaces strictes
		if (sorted_s != sorted_sigma) // O(n) : comparaison de deux vecteurs de taille n
		{
		    F.push_back(sorted_s); // O(n) : copie de sorted_s
		}
		return true; // Continuer la traversée
	};

	// Appel de traverse avec les bons paramètres
	auto tr = st::cofaces<true>(&simplex_tree, simplex_tree.find(sigma));
	traverse(tr, lambda_f);


    return F;
}

// Renvoie le cobord d'un simplexe
// En notant n = sigma.size() et m = cofaces.size() : O(n.log(n) + nm)
vector<simplex_t> MorseSequence::coboundary(const simplex_t& sigma) 
{
	vector<simplex_t> result;
	// Récupération des cofaces
	vector<simplex_t> cofaces_list = this->cofaces(sigma); // O(n.log(n)) : complexité de cofaces

	for (const simplex_t& s : cofaces_list) // O(m) : longueur de cofaces
	{
		if (s.size() == sigma.size() + 1) // O(1)
		{
			result.push_back(s); //O(n) : copie de s
		}
	}

	return result;
}

// Permet de trouver le simplexe v qui est dans B mais pas dans S
// En notant n = S.size() et m = B.size() : O(nm)
std::optional<simplex_t> MorseSequence::find_out(const vector<simplex_t>& B, vector<simplex_t>& S) 
{
	vector<simplex_t> possibilities;

	for (const simplex_t& v : B) // O(m)
	{
		// Vérifie si v n'est pas dans S
		if (std::find(S.begin(), S.end(), v) == S.end()) // O(n)
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
vector<simplex_t> MorseSequence::tri_dim_decroissant(vector<simplex_t>& L) 
{
	std::sort(L.begin(), L.end(), [this](const simplex_t& a, const simplex_t& b) // O(nlog(n)) : complexité de std::sort
	{
		return this->dim(a) > this->dim(b); // Tri décroissant par dimension
	});

	return L;
}

// Pour trier une liste de simplexe par dimensions croissantes
// En notant n = L.size() : O(nlog(n))
vector<simplex_t> MorseSequence::tri_dim_croissant(vector<simplex_t>& L) 
{
	std::sort(L.begin(), L.end(), [this](const simplex_t& a, const simplex_t& b) // O(nlog(n)) : complexité de std::sort
	{
		return this->dim(a) < this->dim(b); // Tri croissant par dimension
	});

	return L;
}


// Renvoie les p-simplexes du complexe simplicial si p est précisé
// Renvoie tous les simplexes si p n'est pas précisé
// En notant n = simplextree.size() : O(n)
vector<simplex_t> MorseSequence::simplices(std::optional<int> p) const 
{
	vector<simplex_t> F;

	// Définir le lambda pour stocker les simplexes
	auto lambda_f = [this, &F](node_ptr, idx_t, const simplex_t& s) -> bool 
	{
		F.push_back(s); // O(1)
		// Continuer la traversée
		return true; 
	};

	// Cas où `p` est None → Traverser tout le complexe
	if (!p) 
	{
		auto tr = st::k_skeleton<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), simplex_tree.dimension()); // O(n) : complexité de st::k_skeleton
		traverse(tr, lambda_f); // O(n) =  O(1) * n : complexité de lambda_f * nombres de simplexes 
	}
	
	// Cas où `p` est défini → Traverser uniquement les `p`-simplexes
	else 
	{
		auto tr = st::k_simplices<true>(&simplex_tree, simplex_tree.find(std::vector<idx_t>{}), *p); // O(m) avec m < n
		traverse(tr, lambda_f);  // O(m) =  O(1) * m : complexité de lambda_f * nombres de simplexes 
	}

	return F;
}



// Création d'une séquence de Morse décroissante
// En notant n = simplextree.size(), m = la taille d'un simplexe : 
std::pair<std::vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>>, int> MorseSequence::morse_seq_dec(const SimplexTree& st)
{
    // Récupération des simplexes triés par dimension décroissante
    vector<simplex_t> K = this->simplices(); // O(n)
    this->tri_dim_decroissant(K); // O(nlog(n))

    // Déclaration et initialisation des variables
    size_t n = K.size();
    size_t i = 0;
    int n_crit = 0;

    std::unordered_map<simplex_t, int> ro;
    vector<simplex_t> L, S;
    vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>> W;

    // Calcul du nombre de cobords pour chaque simplex
    // O(nmlog(m))
    for (const auto& s : K) // O(n)
    { 
        ro[s] = (this->coboundary(s)).size(); // O(mlog(m))
    }

    // Initialisation de L avec les simplexes ayant un seul cobord
    // O(n)
    for (const simplex_t& s : K) // O(n)
    {
        if (ro[s] == 1) // O(1)
        {
            L.push_back(s); // O(1)
        }
    }

    // Boucle principale
    while (i < n) // O(n)
    {
        // Traitement des éléments de L
        while (!L.empty()) // O(n)
        {
            simplex_t tau = L.back(); // O(?)
            L.pop_back(); // O(?)

            std::optional<simplex_t> v_opt = this->find_out(this->coboundary(tau), S); // O(m²log(m))
	    if (!v_opt) // Vérifie si v_opt est nullopt (aucun élément trouvé)
 	    {
    		//std::cout << "Aucun élément trouvé pour tau, on continue.\n";
    		continue;
	    }
	    
	    simplex_t v = *v_opt; // Déballer la valeur

            W.push_back(std::make_pair(v, tau));

            S.push_back(v);
            S.push_back(tau);

            // Construction de la liste des bords
            vector<simplex_t> bords;
            for (const simplex_t& elmt : this->boundary(v))
            {
                bords.push_back(elmt);
            }
            for (const simplex_t& elmt : this->boundary(tau))
            {
                bords.push_back(elmt);
            }

            // Mise à jour de ro et ajout des éléments à L si nécessaire
            for (const simplex_t& mu : bords) 
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
        simplex_t critical = K[i];
        W.push_back(critical);
        n_crit++;

        S.push_back(critical);

        // Mise à jour de ro et ajout des éléments à L si nécessaire
        for (const simplex_t& simplex : this->boundary(critical)) 
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
    //return {W, n_crit};
}


// Création d'une séquence de Morse croissante
// En notant n = simplextree.size() : 
std::pair<std::vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>>, int> MorseSequence::morse_seq_crois(const SimplexTree& K_init)
{
    // Trier les simplexes de K_init par dimension croissante
    vector<simplex_t> K = this->simplices(); // O(n)
    size_t n = K.size();
    size_t i = 0;
    int n_crit = 0;
    
    vector<simplex_t> S, L;
    vector<std::variant<simplex_t, std::pair<simplex_t, simplex_t>>> W;
    std::unordered_map<simplex_t, int> ro;
    
    for (const simplex_t& s : K) {
        ro[s] = 0;
    }
    
    while (i < n - 1) {
        simplex_t sigma = K[i];
        S.push_back(sigma);
        W.push_back(sigma);
        n_crit++;
        
        for (const simplex_t& tau : this->coboundary(sigma)) {
            ro[tau] += 1;
            L.push_back(tau);
        }
        
        this->tri_dim_croissant(L);
        
        while (!L.empty()) {
            simplex_t tau = L.back();
            L.pop_back();
            
            if (ro[tau] == this->dim(tau)) {
                std::optional<simplex_t> v_opt = this->find_out(this->boundary(tau), S);
                if (!v_opt) continue;
                simplex_t v = *v_opt;
                
                W.push_back(std::make_pair(v, tau));
                S.push_back(v);
                S.push_back(tau);
                
                vector<simplex_t> cobords;
                for (const simplex_t& elmt : this->coboundary(v)) {
                    if (std::find(cobords.begin(), cobords.end(), elmt) == cobords.end()) {
                        cobords.push_back(elmt);
                    }
                }
                for (const simplex_t& elmt : this->coboundary(tau)) {
                    if (std::find(cobords.begin(), cobords.end(), elmt) == cobords.end()) {
                        cobords.push_back(elmt);
                    }
                }
                
                for (const simplex_t& mu : cobords) {
                    ro[mu] += 1;
                    if (std::find(L.begin(), L.end(), mu) == L.end()) {
                        L.push_back(mu);
                    }
                }
                
                this->tri_dim_croissant(L);
            }
        }
        
        while (i < n - 1 && std::find(S.begin(), S.end(), K[i]) != S.end()) {
            i++;
        }
    }
    
    return {W, n_crit};
}
