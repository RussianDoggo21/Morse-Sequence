#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/numpy.h>

#include "morse_frame/ref_map.h"
#include "morse_frame/coref_map.h"
#include "morse_sequence/morse_sequence.h"

namespace py = pybind11;

// Conversion py::list -> m_sequence
// Used in _ref_map and _coref_map
m_sequence py_list_to_m_sequence(const py::list& py_W, const SimplexTree& st){
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

py::dict m_frame_to_py_dict(const m_frame& map, const SimplexTree& st) {
    py::dict py_map;

    for (const auto& [key_ptr, bits] : map) {
        // Python key : the full simplex corresponding to key_ptr
        auto vec = st.full_simplex(key_ptr);
        py::tuple py_key(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
            py_key[i] = vec[i];

        // Python value : list of indices where bits are set to 1
        py::list py_val;
        for (std::size_t i = 0; i < bits.size(); ++i) {
            if (bits.test(i))
                py_val.append(py::int_(i));
        }

        py_map[py_key] = py_val;
    }

    return py_map;
}


// Conversion m_sequence -> py::list
// Used in _Min, _Max, _increasing and _decreasing
py::list m_sequence_to_py_list(const m_sequence& W, const SimplexTree& st) {
    std::vector<py::object> result;
    result.reserve(W.size()); // pré-allocation

    for (const auto& v : W) {
        if (std::holds_alternative<node_ptr>(v)) {
            auto sigma = st.full_simplex(std::get<node_ptr>(v));
            result.push_back(py::cast(sigma)); // simplex -> list[int]
        } else {
            auto [a, b] = std::get<std::pair<node_ptr, node_ptr>>(v);
            auto sigma = py::cast(st.full_simplex(a));
            auto tau = py::cast(st.full_simplex(b));
            result.push_back(py::make_tuple(sigma, tau)); // tuple of list[int]
        }
    }

    return py::cast(result); // vector<py::object> → list
}



// Generic function to handle various vector types
template < typename Lambda >
void vector_handler(SimplexTree& st, const py::array_t< idx_t >& simplices, Lambda&& f){
    py::buffer_info s_buffer = simplices.request();

    if (s_buffer.ndim == 1){
        const size_t n = s_buffer.shape[0];
        idx_t* s = static_cast< idx_t* >(s_buffer.ptr);

        for (size_t i = 0; i < n; ++i){
            f(s+i, s+i+1);
        }
    
    } else if (s_buffer.ndim == 2) {

        if (s_buffer.strides[0] <= 0){ return; }
        const size_t d = static_cast< size_t >(s_buffer.shape[1]);
        const size_t n =  static_cast< size_t >(s_buffer.shape[0]);
        idx_t* s = static_cast< idx_t* >(s_buffer.ptr);

        for (size_t i = 0; i < n; ++i){
            f(s+(d*i), s+d*(i+1));
        }
    } else {
        std::cerr << "[vector_handler] ERROR: Unexpected ndim = " << s_buffer.ndim << std::endl;
    }
}


py::tuple _Min_buffered(MorseSequence& ms, const py::array_t<idx_t>& S_buffer, const py::array_t<idx_t>& F_buffer) {
    SimplexTree& st = const_cast<SimplexTree&>(ms.get_simplex_tree());

    node_list cpp_S_full;
    tsl::robin_map<node_ptr, int> cpp_F_full;

    // Fill cpp_S_full with vector_handler
    vector_handler(st, S_buffer, [&](idx_t* b, idx_t* e) {
        simplex_t sigma(b, e);
        node_ptr ptr = st.find(sigma);
        cpp_S_full.push_back(ptr);
    });

    // Fill cpp_F_full with vector_handler (last element = weight)
    vector_handler(st, F_buffer, [&](idx_t* b, idx_t* e) {
        simplex_t sigma(b, e - 1); // everything but the last column 
        int value = *(e - 1);      // last element
        node_ptr ptr = st.find(sigma);
        cpp_F_full[ptr] = value;
    });

    auto [out, ncrit] = ms.Min(cpp_S_full, cpp_F_full);
    return py::make_tuple(m_sequence_to_py_list(out, st), ncrit);
}

py::tuple _Max_buffered(MorseSequence& ms, const py::array_t<idx_t>& S_buffer, const py::array_t<idx_t>& F_buffer) {
    SimplexTree& st = const_cast<SimplexTree&>(ms.get_simplex_tree());

    node_list cpp_S_full;
    tsl::robin_map<node_ptr, int> cpp_F_full;

    // Fill cpp_S_full with vector_handler
    vector_handler(st, S_buffer, [&](idx_t* b, idx_t* e) {
        simplex_t sigma(b, e);
        node_ptr ptr = st.find(sigma);
        cpp_S_full.push_back(ptr);
    });

    // Fill cpp_F_full with vector_handler (last element = weight)
    vector_handler(st, F_buffer, [&](idx_t* b, idx_t* e) {
        simplex_t sigma(b, e - 1); // everything but the last column 
        int value = *(e - 1);      // last element
        node_ptr ptr = st.find(sigma);
        cpp_F_full[ptr] = value;
    });

    auto [out, ncrit] = ms.Max(cpp_S_full, cpp_F_full);
    return py::make_tuple(m_sequence_to_py_list(out, st), ncrit);
}


py::tuple _increasing(MorseSequence& ms, const SimplexTree& st){
    
    // Call of the C++ function
    auto [output, n] = ms.increasing(st);

    // Conversion  m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(output, st);

    return py::make_tuple(out_list, n);
}

py::tuple _decreasing(MorseSequence& ms, const SimplexTree& st){
    
    // Call of the C++ function
    auto [output, n] = ms.decreasing(st);

    // Conversion  m_sequence -> py::list
    py::list out_list = m_sequence_to_py_list(output, st);;

    return py::make_tuple(out_list, n);
}

py::dict _ref_map(MorseSequence& ms, const py::list& py_W){

    // Access to the SimplexTree
    const SimplexTree &st = ms.get_simplex_tree();

    // Conversion py::list -> m_sequence
    m_sequence W = py_list_to_m_sequence(py_W, st);

    // Creation of a Morse Frame
    RefMap ref_map = RefMap(ms, W);
    m_frame bitarray = ref_map.get_bitarray();

    // Conversion m_frame -> py::list
    py::dict py_ref_map = m_frame_to_py_dict(bitarray, st); 

    return py_ref_map;
}

py::dict _coref_map(MorseSequence& ms, const py::list& py_W){

    // Access to the SimplexTree
    const SimplexTree &st = ms.get_simplex_tree();

    // Conversion py::list -> m_sequence
    m_sequence W = py_list_to_m_sequence(py_W, st);

    // Creation of a Morse Frame
    CorefMap ref_map = CorefMap(ms, W);
    m_frame bitarray = ref_map.get_bitarray();

    // Conversion m_frame -> py::list
    py::dict py_coref_map = m_frame_to_py_dict(bitarray, st); 

    return py_coref_map;
}
 

namespace {
    using boundary_fn_1 = node_list (MorseSequence::*)(const node_ptr&) ;
    using boundary_fn_2 = node_list (MorseSequence::*)(const node_ptr&, const tsl::robin_map<node_ptr, bool>&);

    using coboundary_fn_1 = node_list (MorseSequence::*)(const node_ptr&) ;
    using coboundary_fn_2 = node_list (MorseSequence::*)(const node_ptr&, const tsl::robin_map<node_ptr, bool>&);

    using nbboundary_fn_1 = int (MorseSequence::*)(const node_ptr&);
    using nbboundary_fn_2 = int (MorseSequence::*)(const node_ptr&, const tsl::robin_map<node_ptr, bool>&);

    using nbcoboundary_fn_1 = int (MorseSequence::*)(const node_ptr&);
    using nbcoboundary_fn_2 = int (MorseSequence::*)(const node_ptr&, const tsl::robin_map<node_ptr, bool>&);
}



PYBIND11_MODULE(_core, m) {
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
        .def("Max", _Max_buffered)
        .def("Min", _Min_buffered)
        .def("decreasing", _decreasing)
        .def("increasing", _increasing)
        .def("reference_map", _ref_map)
        .def("coreference_map", _coref_map)
        // FONCTION PRINT A DEFINIR DIRECTEMENT EN PYTHON ?? 
        // morse_sequence.py ?
        //.def("print_m_sequence", ) 
        //.def("print_m_frame",)
        ;
}
