#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "morse_sequence.h"
using simplex_t = SimplexTree::simplex_t;
using node_ptr = SimplexTree::node*;
using m_sequence = std::vector<std::variant<node_ptr, std::pair<node_ptr, node_ptr>>>;
using node_list = std::vector<node_ptr>;
using morse_frame = std::unordered_map<node_ptr, node_list>;

namespace py = pybind11;

// Conversion py::list -> node_list
// Used in _Min and _Max
node_list py_list_to_node_list(py::list py_S, const SimplexTree& st){
    node_list S;
    for (auto item : py_S) {
        auto simplex = item.cast<simplex_t>();
        S.push_back(st.find(simplex));  // conversion implicite
    }
    return S;
}

// Conversion py::dict -> std::unordered_map<node_ptr, int>
// Used in _Min and _Max
std::unordered_map<node_ptr, int> py_dict_to_cpp_dict(py::dict py_F, const SimplexTree& st){
    std::unordered_map<node_ptr, int> F;
    for (auto item : py_F) {
        auto key = item.first.cast<simplex_t>();
        auto val = item.second.cast<int>();
        F[st.find(key)] = val;
    }
    return F;
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

// Conversion py::list -> m_sequence
// Used in _ref_map and _coref_map
m_sequence py_list_to_m_sequence(py::list py_W, const SimplexTree& st){
    m_sequence W;      

    for (py::handle item : py_W) {

        /* --- case « critical simplexe » : a single list  ------------------- */
        if (py::isinstance<py::list>(item)) {
            std::vector<int> verts = item.cast<std::vector<int>>();
            node_ptr sigma = st.find(verts);  // Search for the node
            W.emplace_back(sigma);
            continue;
        }

        /* --- cas « free pair (σ,τ) » : a tuple of two lists --------- */
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

// Conversion morse_frame -> py::dict
// Used in _ref_map and _coref_map
py::dict morse_frame_to_py_dict(morse_frame map, const SimplexTree& st){
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

py::tuple _Min(MorseSequence &ms, py::list py_S, py::dict py_F) {
    
    // Access to the SimplexTree
    const SimplexTree& st = ms.get_simplex_tree();
    
    // Conversion py::list -> node_list
    node_list S = py_list_to_node_list(py_S, st);
    
    // Conversion py::dict -> std::unordered_map<node_ptr, int>
    std::unordered_map<node_ptr, int> F = py_dict_to_cpp_dict(py_F, st);

    // Call of the C++ function
    auto [output, n] = ms.Min(S, F);

    // Conversion m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(output, st);

    return py::make_tuple(out_list, n);
}

py::tuple _Max(MorseSequence &ms, py::list py_S, py::dict py_F) {
    
    // Access to the SimplexTree
    const SimplexTree& st = ms.get_simplex_tree();
    
    // Conversion py::list -> node_list
    node_list S = py_list_to_node_list(py_S, st);

    // Conversion py::dict -> std::unordered_map<node_ptr, int>
    std::unordered_map<node_ptr, int> F = py_dict_to_cpp_dict(py_F, st);

    // Call of the C++ function
    auto [output, n] = ms.Max(S, F);

    // Conversion  m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(output, st);

    return py::make_tuple(out_list, n);
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
    morse_frame reference_map = ms.reference_map(W);

    // Conversion morse_frame -> py::list
    py::dict py_ref_map = morse_frame_to_py_dict(reference_map, st);

    return py_ref_map;
}

py::dict _coref_map(MorseSequence &ms, py::list py_W){
    
    // Access to the SimplexTree
    const SimplexTree &st = ms.get_simplex_tree();

    // Conversion py::list -> m_sequence
    m_sequence W = py_list_to_m_sequence(py_W, st);

    // Call of C++ function
    morse_frame coreference_map = ms.coreference_map(W);

    // Conversion morse_frame -> py::list
    py::dict py_coref_map = morse_frame_to_py_dict(coreference_map, st);

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
        //.def("boundary", &MorseSequence::boundary)
        //.def("coboundary", &MorseSequence::coboundary)
        //.def("nbboundary", &MorseSequence::nbboundary)
        //.def("nbcoboundary", &MorseSequence::nbcoboundary)

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
}
