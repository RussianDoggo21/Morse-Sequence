#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "morse_sequence.h"

namespace py = pybind11;

PYBIND11_MODULE(morse_sequence, m) {
    m.doc() = "Interface Python pour MorseSequence";

    py::class_<MorseSequence>(m, "MorseSequence")
        .def(py::init<const SimplexTree&>())
        .def("boundary", &MorseSequence::boundary)
        .def("coboundary", &MorseSequence::coboundary)
        .def("nbboundary", &MorseSequence::nbboundary)
        .def("nbcoboundary", &MorseSequence::nbcoboundary)
        .def("simplices", &MorseSequence::simplices)
        .def("Max", &MorseSequence::Max)
        .def("Min", &MorseSequence::Min)
        .def("find_node", &MorseSequence::find_node);
}
