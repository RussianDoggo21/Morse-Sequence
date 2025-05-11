#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "morse_sequence.h"
using simplex_t = SimplexTree::simplex_t;
using node_ptr = SimplexTree::node*;

namespace py = pybind11;

py::tuple _Min(MorseSequence &ms, py::list py_S, py::dict py_F) {
    
    const SimplexTree& st = ms.get_simplex_tree();
    
    // Conversion py::list -> std::vector<node_ptr>
    std::vector<node_ptr> S;
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
    
    // Conversion py::list -> std::vector<node_ptr>
    std::vector<node_ptr> S;
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

py::tuple _crois(MorseSequence &ms, SimplexTree& st){
    auto [output, n] = ms.morse_seq_crois(st);

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

py::tuple _decrois(MorseSequence &ms, SimplexTree& st){
    auto [output, n] = ms.morse_seq_decrois(st);

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


PYBIND11_MODULE(morse_sequence, m) {
    m.doc() = "Interface Python pour MorseSequence";

    py::class_<MorseSequence>(m, "MorseSequence")
        .def(py::init<const SimplexTree&>())
        .def("boundary", &MorseSequence::boundary)
        .def("coboundary", &MorseSequence::coboundary)
        .def("nbboundary", &MorseSequence::nbboundary)
        .def("nbcoboundary", &MorseSequence::nbcoboundary)
        .def("simplices", &MorseSequence::simplices)
        .def("Max", _Max)
        .def("Min", _Min)
        .def("morse_seq_decrois", _decrois)
        .def("morse_seq_crois", _crois)
        ;
}
