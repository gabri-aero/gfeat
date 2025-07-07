#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fft>
#include <gfeat>
#include <iostream>

namespace py = pybind11;

void init_utils(py::module &m) {
    // Compute repeating orbit
    m.def("sma_from_rgt", &sma_from_rgt);
    // Compute orbital perturbations from J2
    m.def("aop_dot", &aop_dot);
    m.def("ma_dot", &ma_dot);
    m.def("raan_dot", &raan_dot);
    // Rotation matrices
    m.def("R1", &R1);
    m.def("R2", &R2);
    m.def("R3", &R3);
    // Conversion
    m.def("cart2sph", &cart2sph);
    m.def("sph2cart", &sph2cart);
    // ECI to ECRF conversions
    m.def("eci2ecrf",
          static_cast<Eigen::Vector3d (*)(Eigen::Vector3d, double)>(&eci2ecrf),
          py::arg("r_eci"), py::arg("t"));
    m.def("eci2ecrf",
          static_cast<std::pair<Eigen::Vector3d, Eigen::Vector3d> (*)(
              Eigen::Vector3d, Eigen::Vector3d, double)>(&eci2ecrf),
          py::arg("r_eci"), py::arg("v_eci"), py::arg("t"));
    // ECRF to ECI converisons
    m.def("ecrf2eci",
          static_cast<Eigen::Vector3d (*)(Eigen::Vector3d, double)>(&ecrf2eci),
          py::arg("r_ecrf"), py::arg("t"));
    m.def("ecrf2eci",
          static_cast<std::pair<Eigen::Vector3d, Eigen::Vector3d> (*)(
              Eigen::Vector3d, Eigen::Vector3d, double)>(&ecrf2eci),
          py::arg("r_ecrf"), py::arg("v_ecrf"), py::arg("t"));
    // Real FFT
    m.def("rfft", &rfft<double>);
    // Logger
    py::class_<Logger>(m, "Logger")
        .def_static("instance", &Logger::instance,
                    py::return_value_policy::reference)
        .def_property("verbosity", &Logger::get_verbosity,
                      &Logger::set_verbosity);

    py::enum_<Verbosity>(m, "Verbosity")
        .value("Silent", Verbosity::Silent)
        .value("Info", Verbosity::Info);
}
