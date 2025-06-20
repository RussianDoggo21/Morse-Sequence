#include <pybind11/embed.h>   //  ⇐ embeddable interpreter
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>

#include "morse_sequence.h"
#include "../simplex_batch/batch.h"

using simplex_t = SimplexTree::simplex_t;
using node_ptr = SimplexTree::node*;
using m_sequence = std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>;
using node_list = std::vector<node_ptr>;
using m_frame = std::unordered_map<node_ptr, node_list>;
using SimplexList = std::vector<simplex_t>;  // Vector of simplices

namespace py = pybind11;

/*
// Conversion py::array -> node_list
// Used in _Min and _Max
node_list py_list_to_node_list(const py::array& py_S, SimplexTree& st)
{
    node_list S;

    py::array_t<idx_t, py::array::c_style | py::array::forcecast> arr =
        py::cast<py::array>(obj);

    // réserve directement : n = arr.shape(0) 
    S.reserve(arr.request().shape[0]);

    vector_handler(
        st, arr,
        [&](const idx_t* beg, const idx_t* end) {
            S.push_back(st.find(beg, end));       // find(begin,end)
        });
    return S;
}


// Conversion py::dict -> std::unordered_map<node_ptr, int>
// Used in _Min and _Max
std::unordered_map<node_ptr,int>
py_dict_to_cpp_dict(const py::object& obj, SimplexTree& st)
{
    std::unordered_map<node_ptr,int> F;

    // --- numpy array Nx2 -------------------------------------------- 
    if (py::isinstance<py::array>(obj)) {
        auto arr =
          py::array_t<idx_t, py::array::c_style | py::array::forcecast>(obj);
        py::buffer_info info = arr.request();
        if (info.ndim != 2 || info.shape[1] != 2)
            throw std::runtime_error("F array must be Nx2");

        const idx_t* data = static_cast<idx_t*>(info.ptr);
        const std::size_t n = info.shape[0];
        F.reserve(n);

        for (std::size_t i = 0; i < n; ++i) {
            idx_t v   = data[2 * i];
            int   w   = static_cast<int>(data[2 * i + 1]);
            node_ptr p = st.find(&v, &v + 1);
            F.emplace(p, w);
        }
        return F;
    }

    // --- dict Python {simplex : poids} ------------------------------ *
    py::dict d = py::cast<py::dict>(obj);
    F.reserve(d.size());
    for (auto [k, v] : d) {
        simplex_t simplex = k.cast<simplex_t>();
        F.emplace(st.find(simplex), v.cast<int>());
    }
    return F;
}
*/
// Conversion py::list -> m_sequence
// Used in _ref_map and _coref_map
m_sequence py_list_to_m_sequence(py::list py_W, const SimplexTree& st){
    m_sequence W;      

    for (py::handle item : py_W) {

        // --- case « critical simplexe » : a single list  ------------------- 
        if (py::isinstance<py::list>(item)) {
            std::vector<int> verts = item.cast<std::vector<int>>();
            node_ptr sigma = st.find(verts);  // Search for the node
            W.emplace_back(sigma);
            continue;
        }

        // -- case « free pair (σ,τ) » : a tuple of two lists --------- 
        if (py::isinstance<py::tuple>(item) && py::len(item) == 2) {
            py::tuple pair = item.cast<py::tuple>();
            std::vector<int> v1 = pair[0].cast<std::vector<int>>();
            std::vector<int> v2 = pair[1].cast<std::vector<int>>();
            node_ptr sigma = st.find(v1);
            node_ptr tau   = st.find(v2);
            W.emplace_back(std::make_pair(sigma, tau));
            continue;
        }

        throw std::runtime_error("Element of W is neither a simplex "
                                 "nor a (sigma,tau) pair.");
    }
    return W;
}

// Conversion m_frame -> py::dict
// Used in _ref_map and _coref_map
py::dict m_frame_to_py_dict(m_frame map, const SimplexTree& st){
    py::dict py_map;

    for (const auto &[key_ptr, lst] : map) {

        // Python key : the full simplex linked to the node_ptr key_ptr
        auto vec = st.full_simplex(key_ptr);
        py::tuple py_key(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
            py_key[i] = vec[i];
        //py::list py_key = py::cast(py::tuple(py::cast(vec)));

        // Python value : list of critical simplices (or None)
        py::list py_val;
        for (node_ptr v : lst) {
            if (v == nullptr)
                py_val.append(py::none());
            else
                py_val.append(st.full_simplex(v));
        }

        py_map[py_key] = py_val;
    }
    return py_map;
}


// Conversion m_sequence -> py::list
// Used in _Min, _Max, _increasing and _decreasing
py::list m_sequence_to_py_list(m_sequence W, const SimplexTree& st){
    py::list py_W;
    for (const auto& v : W) {
        if (std::holds_alternative<node_ptr>(v)) {
            auto simplex = st.full_simplex(std::get<node_ptr>(v));
            py_W.append(simplex);
        } else {
            auto [a, b] = std::get<std::pair<node_ptr, node_ptr>>(v);
            auto sigma = st.full_simplex(a);
            auto tau = st.full_simplex(b);
            py_W.append(py::make_tuple(sigma, tau));
        }
    }
    return py_W;
}



/*
py::tuple _Min(MorseSequence& ms, py::object py_S, py::object py_F)
{
    // Access to the SimplexTree
    SimplexTree& st = const_cast<SimplexTree&>(ms.get_simplex_tree());

    // Conversion py::list -> node_list
    node_list S = py_list_to_node_list(py_S, st);

    // Conversion py::dict -> std::unordered_map<node_ptr, int>
    std::unordered_map<node_ptr,int> F = py_dict_to_cpp_dict(py_F, st);

    // Call of the C++ function
    auto [out, ncrit] = ms.Min(S, F);

    return py::make_tuple(m_sequence_to_py_list(out, st), ncrit);
}
*/

py::tuple _Min(MorseSequence& ms, const SimplexBatch& batch_S, const SimplexBatch& batch_F){

    // Retrieving of the parameters needed for the C++ function
    std::unordered_map<node_ptr,int> F = batch_to_weight_map(batch_F);
    node_list S = batch_S.nodes;

    // Call of the C++ function
    auto [out, ncrit] = ms.Min(S, F);

    // Conversion of the result into python object
    const SimplexTree& st = ms.get_simplex_tree();
    py::list out_list = m_sequence_to_py_list(out, st);

    return py::make_tuple(out_list, ncrit);
}

/*
py::tuple _Max(MorseSequence& ms, py::object py_S, py::object py_F)
{
    // Access to the SimplexTree
    SimplexTree& st = const_cast<SimplexTree&>(ms.get_simplex_tree());

    // Conversion py::list -> node_list
    node_list S = py_list_to_node_list(py_S, st);

    // Conversion py::dict -> std::unordered_map<node_ptr, int>
    std::unordered_map<node_ptr,int> F = py_dict_to_cpp_map(py_F, st);

    // Call of the C++ function
    auto [out, ncrit] = ms.Max(S, F);

    // Conversion  m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(out, st);

    return py::make_tuple(out_list, ncrit);
}
*/


py::tuple _Max(MorseSequence& ms, const SimplexBatch& batch_S, const SimplexBatch& batch_F){

    // Retrieving of the parameters needed for the C++ function
    std::unordered_map<node_ptr,int> F = batch_to_weight_map(batch_F);
    node_list S = batch_S.nodes;

    // Call of the C++ function
    auto [out, ncrit] = ms.Max(S, F);

    // Conversion of the result into python object
    const SimplexTree& st = ms.get_simplex_tree();
    py::list out_list = m_sequence_to_py_list(out, st);

    return py::make_tuple(out_list, ncrit);
}


py::tuple _increasing(MorseSequence &ms, SimplexTree& st){
    
    // Call of the C++ function
    auto [output, n] = ms.increasing(st);

    // Conversion  m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(output, st);

    return py::make_tuple(out_list, n);
}

py::tuple _decreasing(MorseSequence &ms, SimplexTree& st){
    
    // Call of the C++ function
    auto [output, n] = ms.decreasing(st);

    // Conversion  m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(output, st);;

    return py::make_tuple(out_list, n);
}

py::dict _ref_map(MorseSequence &ms, py::list py_W){
    
    // Access to the SimplexTree
    const SimplexTree &st = ms.get_simplex_tree();

    // Conversion py::list -> m_sequence
    m_sequence W = py_list_to_m_sequence(py_W, st);

    // Call of C++ function
    m_frame reference_map = ms.reference_map(W);

    // Conversion m_frame -> py::list
    py::dict py_ref_map = m_frame_to_py_dict(reference_map, st);

    return py_ref_map;
}

py::dict _coref_map(MorseSequence &ms, py::list py_W){
    
    // Access to the SimplexTree
    const SimplexTree &st = ms.get_simplex_tree();

    // Conversion py::list -> m_sequence
    m_sequence W = py_list_to_m_sequence(py_W, st);

    // Call of C++ function
    m_frame coreference_map = ms.coreference_map(W);

    // Conversion m_frame -> py::list
    py::dict py_coref_map = m_frame_to_py_dict(coreference_map, st);

    return py_coref_map;
}

namespace {
    using boundary_fn_1 = node_list (MorseSequence::*)(node_ptr) ;
    using boundary_fn_2 = node_list (MorseSequence::*)(node_ptr, const std::unordered_map<node_ptr, bool>&);

    using coboundary_fn_1 = node_list (MorseSequence::*)(node_ptr) ;
    using coboundary_fn_2 = node_list (MorseSequence::*)(node_ptr, const std::unordered_map<node_ptr, bool>&);

    using nbboundary_fn_1 = int (MorseSequence::*)(node_ptr);
    using nbboundary_fn_2 = int (MorseSequence::*)(node_ptr, const std::unordered_map<node_ptr, bool>&);

    using nbcoboundary_fn_1 = int (MorseSequence::*)(node_ptr);
    using nbcoboundary_fn_2 = int (MorseSequence::*)(node_ptr, const std::unordered_map<node_ptr, bool>&);
}



PYBIND11_MODULE(morse_sequence, m) {
    m.doc() = "Python interface for MorseSequence";

    py::class_<MorseSequence>(m, "MorseSequence")
        .def(py::init<const SimplexTree&>())

        // Boundary methods
        .def("boundary", static_cast<boundary_fn_1>(&MorseSequence::boundary))
        .def("boundary2", static_cast<boundary_fn_2>(&MorseSequence::boundary))

        // Coboundary methods
        .def("coboundary", static_cast<coboundary_fn_1>(&MorseSequence::coboundary))
        .def("coboundary2", static_cast<coboundary_fn_2>(&MorseSequence::coboundary))

        // Nbboundary methods
        .def("nbboundary", static_cast<nbboundary_fn_1>(&MorseSequence::nbboundary))
        .def("nbboundary2", static_cast<nbboundary_fn_2>(&MorseSequence::nbboundary))

        // Nbcoboundary methods
        .def("nbcoboundary", static_cast<nbcoboundary_fn_1>(&MorseSequence::nbcoboundary))
        .def("nbcoboundary2", static_cast<nbcoboundary_fn_2>(&MorseSequence::nbcoboundary))

        .def("simplices", &MorseSequence::simplices)
        .def("Max", _Max)
        .def("Min", _Min)
        .def("decreasing", _decreasing)
        .def("increasing", _increasing)
        .def("reference_map", _ref_map)
        .def("coreference_map", _coref_map)
        ;

    py::class_<SimplexBatch>(m, "SimplexBatch")
        .def_static("from_python", &SimplexBatch::from_python,
                    py::arg("obj"), py::arg("tree"))
        .def_readonly("simplices", &SimplexBatch::simplices)
        .def_readonly("nodes",     &SimplexBatch::nodes)
        .def_readonly("weights",   &SimplexBatch::weights)
        ;

}

/*
int main() {
    py::scoped_interpreter guard{};          // start / stop Python

    printf("test_bindings_Max: start\n\n");



    SimplexTree st;
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
    for (simplex_t s : L) st.insert(s);

    MorseSequence ms(st);

    //   Build Python objects S (list) and F (dict)

    py::list py_S;
    py::dict py_F;

    //  Example: put every simplex into S and weight 1 into F   
    for ( node_ptr v : ms.simplices(std::nullopt) ) {   
        simplex_t sigma = st.full_simplex(v);   // std::vector<int>

        // 1) convertir en liste Python pour l’appender tel quel 
        py::list py_sigma = py::cast(sigma);
        py_S.append(py_sigma);

        // 2) même contenu, mais immuable pour servir de clé    
        py::tuple key(py_sigma.size());
        for (std::size_t i = 0; i < py_sigma.size(); ++i)
            key[i] = py_sigma[i];

        py_F[key] = 1;           // poids = 1
    }
-
    // Call the pybind11 wrapper _Max
       
    py::tuple result = _Max(ms, py_S, py_F);

    // Unpack Python tuple (m_sequence, int)
    py::list   py_mseq = result[0];
    int        n_crit  = result[1].cast<int>();

    printf("n_crit returned by _Max : %d\n", n_crit);
    printf("m_sequence length      : %zu\n", py_mseq.size());

    printf("\nDone.\n");
    return 0;
}
*/


/*
g++ bindings.cpp morse_sequence.cpp \
    -I/usr/include/python3.13 \
    -I/home/kali/Téléchargements/Morse-Sequence/venv/lib/python3.13/site-packages/pybind11/include \
    -std=c++20 -O3 -Wall -lpython3.13 \
    -o test_bindings
*/