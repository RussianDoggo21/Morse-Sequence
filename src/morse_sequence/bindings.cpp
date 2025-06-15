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

py::tuple _Min(MorseSequence &ms, py::list py_S, py::dict py_F) {
    
    const SimplexTree& st = ms.get_simplex_tree();
    
    // Conversion py::list -> node_list
    node_list S;
    for (auto item : py_S) {
        auto simplex = item.cast<simplex_t>();
        S.push_back(st.find(simplex));  // conversion implicite
    }

    // Conversion py::dict -> std::unordered_map<node_ptr, int>
    std::unordered_map<node_ptr, int> F;
    for (auto item : py_F) {
        auto key = item.first.cast<simplex_t>();
        auto val = item.second.cast<int>();
        F[st.find(key)] = val;
    }

    // Appel de la fonction C++
    auto [output, n] = ms.Min(S, F);

    // Conversion de la sortie -> py::list
    py::list out_list;
    for (const auto& v : output) {
        if (std::holds_alternative<node_ptr>(v)) {
            auto simplex = st.full_simplex(std::get<node_ptr>(v));
            out_list.append(simplex);
        } else {
            auto [a, b] = std::get<std::pair<node_ptr, node_ptr>>(v);
            auto sigma = st.full_simplex(a);
            auto tau = st.full_simplex(b);
            out_list.append(py::make_tuple(sigma, tau));
        }
    }

    return py::make_tuple(out_list, n);
}

py::tuple _Max(MorseSequence &ms, py::list py_S, py::dict py_F) {
    
    const SimplexTree& st = ms.get_simplex_tree();
    
    // Conversion py::list -> node_list
    node_list S;
    for (auto item : py_S) {
        auto simplex = item.cast<simplex_t>();
        S.push_back(st.find(simplex));  // conversion implicite
    }

    // Conversion py::dict -> std::unordered_map<node_ptr, int>
    std::unordered_map<node_ptr, int> F;
    for (auto item : py_F) {
        auto key = item.first.cast<simplex_t>();
        auto val = item.second.cast<int>();
        F[st.find(key)] = val;
    }

    // Appel de la fonction C++
    auto [output, n] = ms.Max(S, F);

    // Conversion de la sortie -> py::list
    py::list out_list;
    for (const auto& v : output) {
        if (std::holds_alternative<node_ptr>(v)) {
            auto simplex = st.full_simplex(std::get<node_ptr>(v));
            out_list.append(simplex);
        } else {
            auto [a, b] = std::get<std::pair<node_ptr, node_ptr>>(v);
            auto sigma = st.full_simplex(a);
            auto tau = st.full_simplex(b);
            out_list.append(py::make_tuple(sigma, tau));
        }
    }

    return py::make_tuple(out_list, n);
}

py::tuple _increasing(MorseSequence &ms, SimplexTree& st){
    auto [output, n] = ms.increasing(st);

    py::list out_list;
    for (const auto& v : output) {
        if (std::holds_alternative<node_ptr>(v)) {
            auto simplex = st.full_simplex(std::get<node_ptr>(v));
            out_list.append(simplex);
        } else {
            auto [a, b] = std::get<std::pair<node_ptr, node_ptr>>(v);
            auto sigma = st.full_simplex(a);
            auto tau = st.full_simplex(b);
            out_list.append(py::make_tuple(sigma, tau));
        }
    }

    return py::make_tuple(out_list, n);
}

py::tuple _decreasing(MorseSequence &ms, SimplexTree& st){
    auto [output, n] = ms.decreasing(st);

    py::list out_list;
    for (const auto& v : output) {
        if (std::holds_alternative<node_ptr>(v)) {
            auto simplex = st.full_simplex(std::get<node_ptr>(v));
            out_list.append(simplex);
        } else {
            auto [a, b] = std::get<std::pair<node_ptr, node_ptr>>(v);
            auto sigma = st.full_simplex(a);
            auto tau = st.full_simplex(b);
            out_list.append(py::make_tuple(sigma, tau));
        }
    }

    return py::make_tuple(out_list, n);
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
    m.doc() = "Interface Python pour MorseSequence";

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
        .def("ms_decreasing", _decreasing)
        .def("ms_increasing", _increasing)
        ;
}
