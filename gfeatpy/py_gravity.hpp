#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gfeat>
#include <iostream>

namespace py = pybind11;

void init_gravity(py::module &m) {
    auto base_functional = py::class_<BaseFunctional>(m, "BaseFunctional",
                                                      R"doc(
            BaseFunctional is an abstract class that serves as parent class to implement 
            different gravity field functionals. For this purpose, it contains a 
            common_degree_factor attribute that is overloaded by the child classes. 

            It is employed internally in both the Global Spherical Harmonics Synthesis 
            and the computation of degree variances. This allows for code modularity
            and simplifies the implementation of the different functionals.
        )doc");

    auto gravity_anomaly =
        py::class_<GravityAnomaly, BaseFunctional>(m, "GravityAnomaly",
                                                   R"doc(
            GravityAnomaly class defines the degree common terms for computation of gravity anomalies
            (in milligals).

            .. math::

                f(l) = 10^5 \frac{\mu}{a_e^2} (l-1)
            )doc");
    auto geoid_height =
        py::class_<GeoidHeight, BaseFunctional>(m, "GeoidHeight",
                                                R"doc(
            GeoidHeight class defines the degree common terms for computation of geoid height.

            .. math::

                f(l) = a_e
            )doc");
    auto equivalent_water_height =
        py::class_<EquivalentWaterHeight, BaseFunctional>(
            m, "EquivalentWaterHeight",
            R"doc(
            EquivalentWaterHeight class defines the degree common terms for computation of EWH 
            (Wahr, 1998).

            .. math::

                f(l) = \frac{a_e \rho_c W_l}{3 \rho_w} \frac{2l+1}{1+k_l'}

            The class also include the possibility to apply Gaussian smoothing with an
            input smoothing radius (Jekeli, 1981). The computation of the SH coefficients of 
            the averaging function :math:`W_l` follows the continuous fraction approach proposed by 
            Piretzidis (2019).

            Inherits from:
                BaseFunctional        
        )doc");

    gravity_anomaly.def(py::init<>(),
                        R"doc(
            Constructor for GravityAnomaly.
                
            Example
            --------
            Create a BaseFunctional instance for gravity anomalies

                gravity_anomaly = GravityAnomaly()
    )doc");
    geoid_height.def(py::init<>(),
                     R"doc(
            Constructor for GeoidHeight.
                
            Example
            --------
            Create a BaseFunctional instance for geoid heights

                geoid_height = GeoidHeight()
    )doc");
    equivalent_water_height.def(py::init<double>(), py::arg("smoothing_radius"),
                                R"doc(
            Constructor for EquivalentWaterHeight.

            Parameters
            -----------
            smoothing_radius : float
                Value at which the Gaussian spatial smoothing kernel decays to 1/2 of the 
                initial value.
                
            Example
            --------
            Create an instance with 200 km smoothing radius::

                ewh = EquivalentWaterHeight(200e3)
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
