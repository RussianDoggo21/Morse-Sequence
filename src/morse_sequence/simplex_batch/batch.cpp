#include "batch.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using idx_t = std::size_t;


// Copied from _simplextree.cpp in the SimplexTree library
template < typename Lambda >
void vector_handler(SimplexTree& st,const py::array_t<idx_t, py::array::c_style | py::array::forcecast>& arr,Lambda&& f){
    py::buffer_info info = arr.request();
    const idx_t* data    = static_cast<idx_t*>(info.ptr);

    if (info.ndim == 1) {
        const idx_t n = info.shape[0];
        for (idx_t i = 0; i < n; ++i)
            f(data + i, data + i + 1);                // [v,v+1)
    }
    else if (info.ndim == 2) {
        const idx_t n = info.shape[0];
        const idx_t d = info.shape[1];
        for (idx_t i = 0; i < n; ++i)
            f(data + i * d, data + (i + 1) * d);      // [row,row+d)
    }
    else
        throw std::runtime_error("array ndim must be 1 or 2");
}



/*------------------------------------------------------------------*/
SimplexBatch SimplexBatch::from_python(const py::object& obj, SimplexTree& st){

    SimplexBatch batch;

    /* ---------- 1. ndarray fast paths ----------------------------- */
    if (py::isinstance<py::array>(obj))
    {
        auto arr = py::array_t<idx_t,
                   py::array::c_style | py::array::forcecast>(obj);
        auto inf = arr.request();

        /* ----- case (N,2)  →  F‑sequence [vertex, weight] ------- */
        if (inf.ndim == 2 && inf.shape[1] == 2)
        {
            const idx_t* data = static_cast<idx_t*>(inf.ptr);
            const idx_t n = inf.shape[0];

            batch.simplices.reserve(n);
            batch.nodes.reserve(n);
            batch.weights.reserve(n);

            for (idx_t i = 0; i < n; ++i)
            {
                idx_t v = data[2*i];
                int   w = static_cast<int>(data[2*i+1]);

                simplex_t sig{v};                       // 0‑simplex
                batch.simplices.push_back(sig);
                batch.nodes.push_back(st.find(simplex_t(&v, &v+1)));
                batch.weights.push_back(w);
            }
            return batch;
        }

        /* ----- case 1‑D or 2‑D (N,d)  →  cosimplicial list S ---- */
        vector_handler(st, arr,
            [&](const idx_t* b, const idx_t* e) {
                batch.simplices.emplace_back(b, e);
                batch.nodes.push_back(st.find(simplex_t(b, e)));
            });
        return batch;
    }

    /* ---------- 2. dict { simplex : weight }  -------------------- */
    if (py::isinstance<py::dict>(obj))
    {
        py::dict d = py::cast<py::dict>(obj);
        batch.simplices.reserve(d.size());
        batch.nodes.reserve(d.size());
        batch.weights.reserve(d.size());

        for (auto [k, v] : d)
        {
            simplex_t sig = k.cast<simplex_t>();
            int       w   = v.cast<int>();

            batch.simplices.push_back(sig);
            node_ptr np = st.find(sig);
            batch.nodes.push_back(np);
            batch.weights.push_back(w);
        }
        return batch;
    }

    /* ---------- 3. generic Python list / tuple  ------------------ */
    py::sequence seq = py::cast<py::sequence>(obj);
    batch.simplices.reserve(seq.size());
    batch.nodes.reserve(seq.size());

    for (py::handle h : seq)
    {
        /* list of simplexes  -> (simplex)            */
        if (!py::isinstance<py::tuple>(h) || py::len(h) != 2) {
            simplex_t sig = h.cast<simplex_t>();
            batch.simplices.push_back(sig);
            batch.nodes.push_back(st.find(sig));
            continue;
        }

        /* list of pairs (simplex, weight) ----------- */
        py::tuple tup = py::cast<py::tuple>(h);
        simplex_t sig = tup[0].cast<simplex_t>();
        int       w   = tup[1].cast<int>();

        batch.simplices.push_back(sig);
        batch.nodes.push_back(st.find(sig));
        batch.weights.push_back(w);
    }
    return batch;
}
