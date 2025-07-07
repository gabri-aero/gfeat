#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gfeat>
#include <iostream>

namespace py = pybind11;

void init_gravity(py::module &m) {
    auto base_functional = py::class_<BaseFunctional>(m, "BaseFunctional");

    auto gravity_anomaly =
        py::class_<GravityAnomaly, BaseFunctional>(m, "GravityAnomaly");
    auto geoid_height =
        py::class_<GeoidHeight, BaseFunctional>(m, "GeoidHeight");
    auto equivalent_water_height =
        py::class_<EquivalentWaterHeight, BaseFunctional>(
            m, "EquivalentWaterHeight",
            R"doc(
            EquivalentWaterHeight class defines the common terms for computation of EWH.

            The class also includes the possibility to apply Gaussian smoothing with an
            input smoothing radius.

            Inherits from:
                BaseFunctional

            Example usage:
                ewh = EquivalentWaterHeight(200e3)            
        )doc");

    gravity_anomaly.def(py::init<>());
    geoid_height.def(py::init<>());
    equivalent_water_height.def(py::init<double>(), R"doc(
            Constructor for EquivalentWaterHeight.

            Parameters
            -----------
            smoothing_radius : float
                Value at which the Gaussian spatial smoothing kernel decays to 1/2 of the 
                initial value.
    )doc");

    auto sh = py::class_<SphericalHarmonics>(m, "SphericalHarmonics");

    sh.def(py::init<double>())
        .def("potential", &SphericalHarmonics::potential, py::arg("r_ecrf"))
        .def("gravity_anomaly", &SphericalHarmonics::gravity_anomaly,
             py::arg("r_ecrf"))
        .def("geoid_height", &SphericalHarmonics::geoid_height,
             py::arg("r_ecrf"))
        .def("gravity", &SphericalHarmonics::gravity, py::arg("r_ecrf"))
        .def("synthesis", &SphericalHarmonics::synthesis)
        .def("degree_variance", &SphericalHarmonics::degree_variance)
        .def("rms_per_coefficient_per_degree",
             &SphericalHarmonics::rms_per_coefficient_per_degree)
        .def_readwrite("mu", &SphericalHarmonics::mu)
        .def_readwrite("ae", &SphericalHarmonics::ae)
        .def_readwrite("coefficients", &SphericalHarmonics::coefficients);
    ;
    auto sh_error =
        py::class_<SphericalHarmonicsError>(m, "SphericalHarmonicsError");

    sh_error.def(py::init<int>(), py::arg("l_max"))
        .def(py::init<int, std::string, std::string>(), py::arg("l_max"),
             py::arg("file"), py::arg("root") = DATA_DIR)
        .def("synthesis", &SphericalHarmonicsError::synthesis)
        .def("get_Pxx", &SphericalHarmonicsError::get_Pxx)
        .def("degree_variance", &SphericalHarmonicsError::degree_variance)
        .def("rms_per_coefficient_per_degree",
             &SphericalHarmonicsError::rms_per_coefficient_per_degree)
        .def_readwrite("Pxx", &SphericalHarmonicsError::Pxx);
    ;

    auto gravity_field =
        py::class_<GravityField, SphericalHarmonics>(m, "GravityField");

    gravity_field.def(py::init<std::string, int, int>());

    auto datetime = py::class_<DateTime>(m, "DateTime");

    datetime.def(py::init<>())
        .def(py::init<int, int, int, int, int, int>(), py::arg("year"),
             py::arg("month"), py::arg("day"), py::arg("h") = 0,
             py::arg("m") = 0, py::arg("s") = 0)
        .def_readwrite("year", &DateTime::year)
        .def_readwrite("month", &DateTime::month)
        .def_readwrite("day", &DateTime::day)
        .def_readwrite("h", &DateTime::h)
        .def_readwrite("m", &DateTime::m)
        .def_readwrite("s", &DateTime::s);

    auto aod1b_type = py::enum_<AOD1BType>(m, "AOD1BType");

    aod1b_type.value("ATM", AOD1BType::ATM)
        .value("OCN", AOD1BType::OCN)
        .value("GLO", AOD1BType::GLO)
        .value("OBA", AOD1BType::OBA)
        .export_values();

    auto aod1b = py::class_<AOD1B>(m, "AOD1B");

    aod1b.def(py::init<>())
        .def("load", &AOD1B::load, py::arg("filename"),
             py::arg("root") = DATA_DIR)
        .def("get", &AOD1B::get);
}
