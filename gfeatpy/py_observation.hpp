#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>

#include <gfeat>

namespace py = pybind11;

void init_observation(py::module &m) {

    py::enum_<LongitudePolicy>(m, "LongitudePolicy")
        .value("INTERLEAVING", LongitudePolicy::INTERLEAVING)
        .value("OVERLAPPING", LongitudePolicy::OVERLAPPING);

    auto abstrack_kite_system =
        py::class_<AbstractKiteSystem, std::shared_ptr<AbstractKiteSystem>>(
            m, "AbstractKiteSystem");

    abstrack_kite_system.def("solve", &AbstractKiteSystem::block_solve)
        .def("get_sigma_x", &AbstractKiteSystem::get_sigma_x)
        .def("degree_variance", &AbstractKiteSystem::get_degree_variance)
        .def("synthesis", &AbstractKiteSystem::synthesis)
        .def("synthesis_average", &AbstractKiteSystem::synthesis_average)
        .def("get_N_blocks", &AbstractKiteSystem::get_N_blocks)
        .def("get_H_blocks", &AbstractKiteSystem::get_H_blocks)
        .def("get_N", &AbstractKiteSystem::get_N)
        .def("set_kaula_regularization",
             &AbstractKiteSystem::set_kaula_regularization)
        .def("set_solution_time_window",
             &AbstractKiteSystem::set_solution_time_window);

    auto base_observation =
        py::class_<BaseObservation, AbstractKiteSystem,
                   std::shared_ptr<BaseObservation>>(m, "BaseObservation");

    base_observation
        .def("set_observation_error", &BaseObservation::set_observation_error)
        .def("get_Pyy", &BaseObservation::get_Pyy)
        .def("get_H", &BaseObservation::get_H)
        .def("simulate_observations", &BaseObservation::simulate_observations)
        .def("get_radius", &BaseObservation::get_radius);

    auto line_potential =
        py::class_<Potential, BaseObservation, std::shared_ptr<Potential>>(
            m, "Potential");

    line_potential.def(py::init<int, double, int, int, double, double>(),
                       py::arg("l_max"), py::arg("I"), py::arg("Nr"),
                       py::arg("Nd"), py::arg("we_0") = 0, py::arg("wo_0") = 0);

    auto radial = py::class_<Radial, BaseObservation, std::shared_ptr<Radial>>(
        m, "Radial");

    radial.def(py::init<int, double, int, int, double, double>(),
               py::arg("l_max"), py::arg("I"), py::arg("Nr"), py::arg("Nd"),
               py::arg("we_0") = 0, py::arg("wo_0") = 0);

    auto along_track =
        py::class_<AlongTrack, BaseObservation, std::shared_ptr<AlongTrack>>(
            m, "AlongTrack");

    along_track.def(py::init<int, double, int, int, double, double>(),
                    py::arg("l_max"), py::arg("I"), py::arg("Nr"),
                    py::arg("Nd"), py::arg("we_0") = 0, py::arg("wo_0") = 0);

    auto cross_track =
        py::class_<CrossTrack, BaseObservation, std::shared_ptr<CrossTrack>>(
            m, "CrossTrack");

    cross_track.def(py::init<int, double, int, int, double, double>(),
                    py::arg("l_max"), py::arg("I"), py::arg("Nr"),
                    py::arg("Nd"), py::arg("we_0") = 0, py::arg("wo_0") = 0);

    auto collinear =
        py::class_<Collinear, BaseObservation, std::shared_ptr<Collinear>>(
            m, "Collinear");

    collinear
        .def(py::init<int, double, double, int, int, double, double>(),
             py::arg("l_max"), py::arg("I"), py::arg("eta"), py::arg("Nr"),
             py::arg("Nd"), py::arg("we_0") = 0, py::arg("wo_0") = 0)
        .def("set_observation_error",
             [](Collinear &self, std::function<double(double)> range_asd,
                std::function<double(double)> accelerometer_asd) {
                 self.set_observation_error(range_asd, accelerometer_asd);
             });

    auto multi_observation =
        py::class_<MultiObservation, AbstractKiteSystem,
                   std::shared_ptr<MultiObservation>>(m, "MultiObservation");

    multi_observation
        .def(py::init<int, int, int,
                      std::vector<std::shared_ptr<BaseObservation>>>(),
             py::arg("l_max"), py::arg("Nr"), py::arg("Nd"),
             py::arg("observations"))
        .def("update", &MultiObservation::update)
        .def("get_observations", &MultiObservation::get_observations);

    auto gps =
        py::class_<GPS, MultiObservation, std::shared_ptr<GPS>>(m, "GPS");

    gps.def(py::init<int, double, int, int, double, double>(), py::arg("l_max"),
            py::arg("I"), py::arg("Nr"), py::arg("Nd"), py::arg("we_0"),
            py::arg("wo_0"));

    auto constellation =
        py::class_<Constellation, MultiObservation,
                   std::shared_ptr<Constellation>>(m, "Constellation");

    constellation
        .def(
            py::init<int, Eigen::VectorXd, double, int, int, LongitudePolicy>(),
            py::arg("l_max"), py::arg("I"), py::arg("eta"), py::arg("Nr"),
            py::arg("Nd"),
            py::arg("LongitudePolicy") = LongitudePolicy::INTERLEAVING)
        .def("set_observation_error", &Constellation::set_observation_error);

    auto bender =
        py::class_<Bender, Constellation, std::shared_ptr<Bender>>(m, "Bender");

    bender.def(
        py::init<int, Eigen::VectorXd, double, int, int, LongitudePolicy>(),
        py::arg("l_max"), py::arg("I"), py::arg("eta"), py::arg("Nr"),
        py::arg("Nd"),
        py::arg("LongitudePolicy") = LongitudePolicy::INTERLEAVING);
}
