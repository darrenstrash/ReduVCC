#include <pybind11/pybind11.h>

#include "redu_vcc/redu_vcc.h"
#include "redu_vcc/reducer.h"
#include "branch_and_reduce/b_and_r.h"

namespace py = pybind11;

PYBIND11_MODULE(vcc, m)
{
    py::class_<redu_vcc>(m, "redu_cc")
      .def(py::init<>())
      ;

}