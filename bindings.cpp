#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "morse_sequence.h"

namespace py = pybind11;

PYBIND11_MODULE(morse_sequence, m) {
    m.doc() = "Interface Python pour MorseSequence";

    py::class_<MorseSequence>(m, "MorseSequence")
        .def(py::init<const SimplexTree&>())  // Constructeur avec un SimplexTree
        .def("dim", &MorseSequence::dim)
        .def("boundary", &MorseSequence::boundary)
        .def("coboundary", &MorseSequence::coboundary)
        .def("find_out", &MorseSequence::find_out)
        .def("tri_dim_decroissant", &MorseSequence::tri_dim_decroissant)
        .def("tri_dim_croissant", &MorseSequence::tri_dim_croissant)
        .def("simplices", &MorseSequence::simplices)
        .def("morse_seq_dec", &MorseSequence::morse_seq_dec)
        .def("morse_seq_crois", &MorseSequence::morse_seq_crois);
}
