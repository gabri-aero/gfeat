#include <pybind11/pybind11.h>

#include <gfeat>

#include "py_gravity.hpp"
#include "py_observation.hpp"
#include "py_utils.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "Main module for the GFEAT package";

    py::module observation =
        m.def_submodule("observation", "Observation submodule");
    init_observation(observation);

    py::module gravity = m.def_submodule("gravity", "Gravity submodule");
    init_gravity(gravity);

    py::module utils = m.def_submodule("utils", "Utils submodule");
    init_utils(utils);

    // Define the planet class
    py::class_<Planet>(m, "Planet",
                       "This class defines a global variable that stores the "
                       "properties of the central body that is under analysis")
        .def(py::init<>()) // Default constructor
        .def_readwrite("mu", &Planet::mu,
                       "Standard gravitational parameter [m³/s²]")
        .def_readwrite("theta_dot", &Planet::theta_dot,
                       "Angular rotational velocity [rad/s]")
        .def_readwrite("C20", &Planet::C20,
                       "Unnormalized flattening coefficient")
        .def_readwrite("ae", &Planet::ae,
                       "Equatorial radius, or Brillouin sphere [m]")
        .def_readwrite("rho_e", &Planet::rho_e, "Mean crustal density [kg/m³]")
        .def_readwrite("rho_w", &Planet::rho_w,
                       "Water density at Sea Level for EWH [kg/m³]");

    // Bind the global planet instance
    m.attr("planet") = &planet;
}
