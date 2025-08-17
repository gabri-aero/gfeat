#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>

#include <gfeat>

namespace py = pybind11;

void init_observation(py::module &m) {

    py::enum_<LongitudePolicy>(m, "LongitudePolicy", R"doc(
            Enumeration that defines the two built-in longitude separation policies for the :attr:`~gfeatpy.observation.Constellation`
            class.
        )doc")
        .value("INTERLEAVING", LongitudePolicy::INTERLEAVING, R"doc(
            The ground-tracks are perfectly interleaved at the equator.
        )doc")
        .value("OVERLAPPING", LongitudePolicy::OVERLAPPING, R"doc(
            The ground-tracks overlap at the equator.
        )doc");

    auto abstrack_kite_system =
        py::class_<AbstractKiteSystem, std::shared_ptr<AbstractKiteSystem>>(
            m, "AbstractKiteSystem", R"doc(
            AbstractKiteSystem class defines the base object that handles the block-kite structure 
            of the linear system of equations along a circular repeating ground-track.
        )doc");

    abstrack_kite_system
        .def("solve", &AbstractKiteSystem::block_solve,
             R"doc(
            This function solves for the parameter covariance :math:`P_{xx}` making use of least squares
            error propagation.

            .. math::

                P_{xx} = (H^T P_{yy}^{-1} H + P_{cc}^{-1})^{-1}

            Cholesky decomposition allows for leveraging the symmetry of the normal matrix for efficient
            inversion.
    )doc")
        .def("get_sigma_x", &AbstractKiteSystem::get_sigma_x,
             R"doc(
            Returns the compact matrix storing the standard deviation of the parameters, the SH coefficients. See 
            :attr:`~gfeatpy.gravity.SphericalHarmonics.coefficients` for the matrix structure.
    )doc")
        .def("degree_variance", &AbstractKiteSystem::degree_variance,
             R"doc(
            This function computes degree variance spectrum from the standard deviation of the SH coefficients.

            .. math::
                \sigma_l^2 = \sum_{m=0}^{l} \sigma^2(C_{lm}) + \sigma^2(S_{lm})
    )doc")
        .def("rms_per_coefficient_per_degree",
             &AbstractKiteSystem::rms_per_coefficient_per_degree,
             R"doc(
            This function computes RMS per coefficient per degree from the standard deviation 
            of the SH coefficients.

            .. math::
                \delta_l = \sqrt{\sum_{m=0}^{l} \frac{\sigma^2(C_{lm}) + \sigma^2(S_{lm})}{2l+1}}
            )doc")
        .def("synthesis", &AbstractKiteSystem::synthesis, py::arg("n_lon"),
             py::arg("n_lat"), py::arg("functional"), R"doc(
            This function performs covariance propagation from the parameter covariance into the sphere according to the 
            input gravity field functional.
            
            .. math::

                \sigma_{y}(\lambda, \phi) = v P_{xx} v^T

            It makes use of the two-step fold proposed by Rizos (1979) for efficient computation.
            It does not use FFT since the sparsity of the design and parameter covariance matrices makes it 
            less efficient than summation along the dense blocks.

            Parameters
            -----------
            n_lon : int
                Number of nodes along the longitude direction.
            n_lat : int
                Number of nodes along the latitude direction.
            functional : BaseFunctional
                Gravity field functional to perform the synthesis.
                
            Return
            --------
            numpy.ndarray
                A 2D array with longitude values for the lon/lat grid.
            numpy.ndarray
                A 2D array with latitude values for lon/lat grid.
            numpy.ndarray
                A 2D array with the commission error values for
                the input functional on the lon/lat grid.
             )doc")
        .def("synthesis_average", &AbstractKiteSystem::synthesis_average,
             py::arg("functional"), R"doc(
            This function computes the average on the sphere of the comission error of the input 
            gravity field functional. For this purpose, orthogonal properties of spherical harmonics 
            are leveraged.
            
            .. math::

                \hat{\sigma}_{y} = \frac{1}{4\pi} \int_{\Omega} \sigma_y(\lambda, \phi) \, d\Omega
        
            Parameters
            -----------
            functional : BaseFunctional
                Gravity field functional to perform the synthesis.
                
            Return
            --------
            float
                Commission error spherical average for input functional.

             )doc")
        .def("get_N", &AbstractKiteSystem::get_N,
             R"doc(
            Getter for full normal matrix.

            Return
            --------
            ndarray
                Full normal matrix.
            )doc")
        .def("get_Pxx", &AbstractKiteSystem::get_Pxx,
             R"doc(
            Getter for full parameter covariance matrix.

            Return
            --------
            ndarray
                Full parameter covariance matrix.
            )doc")
        .def("set_kaula_regularization",
             &AbstractKiteSystem::set_kaula_regularization, py::arg("K") = 1e-5,
             R"doc(
            This function introduces Kaula's rule to solve the linear system. It provides a constraint on the 
            standard deviation of the parameters as defined by Kaula's power rule:

            .. math::        
                \sigma\{\bar{C}_{lm}, \bar{S}_{lm}\} = \frac{K}{l^2}     

                
            Parameters
            -----------
            K : float
                Kaula's constant for the empirical Kaula's rule"

            )doc")
        .def("set_solution_time_window",
             &AbstractKiteSystem::set_solution_time_window,
             py::arg("time_window"),
             R"doc(
            Setter for the time window :math:`T` of data accumulation that contributes to a gravity field solution.
            It is therefore the sampling time of the gravity field. It defines the frequency resolution :math:`\Delta f=1/T`
            that is employed to compute the error in the lumped coefficients from the amplitude spectral density.

            By default, the solution time window  is set to the repeatability period of the ground-track. 
            
            Parameters
            -----------
            time_window : float
                Solution time window in days.

                
            .. warning::

                Note that when there is not a commensurability with the ground-track repeatability time, the solutions might
                be inaccurate, since the periodicity of the orbit is not fulfilled. Solutions with a time window lower than 
                the ground-track repeatability time are highly unrecommended.

            )doc");

    auto base_observation = py::class_<BaseObservation, AbstractKiteSystem,
                                       std::shared_ptr<BaseObservation>>(
        m, "BaseObservation", R"doc(
            BaseObservation class defines the base object for any observation along a circular orbit.
            It is a derived class from :attr:`~gfeatpy.observation.AbstractKiteSystem`, which provides the structure of
            the linear system of equations along a circular repeating ground-track. This class defines some additional
            settings that depend on the observations themselves, such as the vertical length of the block-design matrices, 
            the actual orbital radius or the amplitude spectral density for the lumped coefficients. 
            Moreover, it provides a framework to define any arbitrary design matrix.
        )doc");

    base_observation
        .def("set_observation_error", &BaseObservation::set_observation_error,
             py::arg("asd"),
             R"doc(
            
            This function sets the lumped coefficients error in the frequency domain from the Amplitude Spectral Density of the 
            observation (Sneeuw, 2000).

            .. math::

                \sigma_{km} = \int_{f_{km}+\Delta f/2}^{f_{km}+\Delta f/2} A^2(f) \, df \approx A^2(f_{km}) \Delta f

            Parameters
            -----------
            asd : Callable[[float], float]
                Function that defines the amplitude spectral density for an input frequency.
            

            )doc")
        .def("get_Pyy", &BaseObservation::get_Pyy, R"doc(
            
            Getter for the full observation covariance matrix.
            
            Return
            --------
            ndarray
                Full observation covariance matrix.
            
            )doc")
        .def("get_H", &BaseObservation::get_H, R"doc(
            
            Getter for the full design matrix.
            
            Return
            --------
            ndarray
                Full design matrix.
            )doc")
        .def("simulate_observations", &BaseObservation::simulate_observations,
             py::arg("gravity_field"),
             R"doc(
            
            This function simulates analytically the values of the lumped coefficients :math:`A_{km}, B_{km}`
            from an input :attr:`~gfeatpy.gravity.GravityField`. 

            Parameters
            -----------
            gravity_field : GravityField
                Gravity field object to fill the parameter vector :math:`x`.
            
            Return
            --------
            dict[int, tuple[float, float]]
                Dictionary that contains the integered frequencies as keys and the associated :math:`A_{km}, B_{km}`
                coefficients as values.            
            )doc")
        .def("get_radius", &BaseObservation::get_radius, R"doc(
            
            Getter for orbital radius associated to the Repeating Ground-Track selected.
            
            Return
            --------
            float
                Orbital radius associated to the Repeating Ground-Track of the observation.

            )doc");

    auto line_potential =
        py::class_<Potential, BaseObservation, std::shared_ptr<Potential>>(
            m, "Potential", R"doc(
        
        This class defines the perturbing potential observation for a circular orbit.
            
        )doc");

    line_potential.def(py::init<int, int, int, double, double, double>(),
                       py::arg("l_max"), py::arg("Nr"), py::arg("Nd"),
                       py::arg("I"), py::arg("we_0") = 0, py::arg("wo_0") = 0,
                       R"doc(
        
        Constructor for the Radial class.

        Parameters
        -----------
        l_max : int
            Cut-off degree.
        Nr : int
            Number of orbital revolutions per ground-track repeatability period. 
        Nd : int
            Nodal days per ground-track repeatability period. 
        I : float
            Inclination [rad].
        we_0 : float
            Initial Earth-fixed longitude of the ascending node [rad].
        wo_0 : float
            Initial argument of latitude [rad].
        )doc");

    auto radial = py::class_<Radial, BaseObservation, std::shared_ptr<Radial>>(
        m, "Radial", R"doc(
        
        This class defines the position displacement observation in the radial
        direction.
            
        )doc");

    radial.def(py::init<int, int, int, double, double, double>(),
               py::arg("l_max"), py::arg("Nr"), py::arg("Nd"), py::arg("I"),
               py::arg("we_0") = 0, py::arg("wo_0") = 0,
               R"doc(
        
        Constructor for the Radial class.

        Parameters
        -----------
        l_max : int
            Cut-off degree.
        Nr : int
            Number of orbital revolutions per ground-track repeatability period. 
        Nd : int
            Nodal days per ground-track repeatability period. 
        I : float
            Inclination [rad].
        we_0 : float
            Initial Earth-fixed longitude of the ascending node [rad].
        wo_0 : float
            Initial argument of latitude [rad].
        )doc");

    auto along_track =
        py::class_<AlongTrack, BaseObservation, std::shared_ptr<AlongTrack>>(
            m, "AlongTrack", R"doc(
        
        This class defines the position displacement observation in the along-track
        direction.
            
        )doc");

    along_track.def(py::init<int, int, int, double, double, double>(),
                    py::arg("l_max"), py::arg("Nr"), py::arg("Nd"),
                    py::arg("I"), py::arg("we_0") = 0, py::arg("wo_0") = 0,
                    R"doc(
        
        Constructor for the AlongTrack class.

        Parameters
        -----------
        l_max : int
            Cut-off degree.
        Nr : int
            Number of orbital revolutions per ground-track repeatability period. 
        Nd : int
            Nodal days per ground-track repeatability period. 
        I : float
            Inclination [rad].
        we_0 : float
            Initial Earth-fixed longitude of the ascending node [rad].
        wo_0 : float
            Initial argument of latitude [rad].
        )doc");

    auto cross_track =
        py::class_<CrossTrack, BaseObservation, std::shared_ptr<CrossTrack>>(
            m, "CrossTrack", R"doc(
        
        This class defines the position displacement observation in the cross-track
        direction.
            
        )doc");

    cross_track.def(py::init<int, int, int, double, double, double>(),
                    py::arg("l_max"), py::arg("Nr"), py::arg("Nd"),
                    py::arg("I"), py::arg("we_0") = 0, py::arg("wo_0") = 0,
                    R"doc(
        
        Constructor for the CrossTrack class.

        Parameters
        -----------
        l_max : int
            Cut-off degree.
        Nr : int
            Number of orbital revolutions per ground-track repeatability period. 
        Nd : int
            Nodal days per ground-track repeatability period. 
        I : float
            Inclination [rad].
        we_0 : float
            Initial Earth-fixed longitude of the ascending node [rad].
        wo_0 : float
            Initial argument of latitude [rad].
        )doc");

    auto collinear =
        py::class_<Collinear, BaseObservation, std::shared_ptr<Collinear>>(
            m, "Range");

    collinear
        .def(py::init<int, int, int, double, double, double, double>(),
             py::arg("l_max"), py::arg("Nr"), py::arg("Nd"), py::arg("I"),
             py::arg("rho_0"), py::arg("we_0") = 0, py::arg("wo_0") = 0,
             R"doc(
        
        Constructor for the Range class.

        Parameters
        -----------
        l_max : int
            Cut-off degree.
        Nr : int
            Number of orbital revolutions per ground-track repeatability period. 
        Nd : int
            Nodal days per ground-track repeatability period. 
        I : float
            Inclination [rad].
        rho_0 : float
            Nominal intersatellite distance [m].
        we_0 : float
            Initial Earth-fixed longitude of the ascending node [rad].
        wo_0 : float
            Initial argument of latitude [rad].
        )doc")
        .def(
            "set_observation_error",
            [](Collinear &self, std::function<double(double)> range_asd,
               std::function<double(double)> accelerometer_asd) {
                self.set_observation_error(range_asd, accelerometer_asd);
            },
            py::arg("range_instrument_asd"), py::arg("accelerometer_asd"),
            R"doc(
            
            An important error contribution to the range observation error is also
            the error in the accelerometer. This error propagates to the range observation.
            This function accounts for this propagation internally.

            
            Parameters
            -----------
            range_instrument_asd : Callable[[float], float]
                Ranging instrument amplitude spectral density for an input frequency.
            accelerometer_asd : Callable[[float], float]
                Accelerometer amplitude spectral density for an input frequency.

            )doc");

    auto multi_observation = py::class_<MultiObservation, AbstractKiteSystem,
                                        std::shared_ptr<MultiObservation>>(
        m, "MultiObservation", R"doc(
        
        This class takes a list of attr:`~gfeatpy.observation.BaseObservation` objects to compute 
        the combined performance of multiple observations.

        .. warning::
            To maintain the simplicity of the block-kite system, this class only supports observations
            with the same associated ground-track repeatability conditions (same number of revolutions :math:`N_r`
            as well as nodal days :math:`N_d`). In this way, the individual block-kite structure of every observation
            can be applied to the system of combined observations.
        )doc");

    multi_observation
        .def(py::init<int, int, int,
                      std::vector<std::shared_ptr<BaseObservation>>>(),
             py::arg("l_max"), py::arg("Nr"), py::arg("Nd"),
             py::arg("observations"), R"doc(
             
            Constructor for MultiObservation class        

            Parameters
            -----------
            l_max : int
                Cut-off degree.
            Nr : int
                Number of orbital revolutions per ground-track repeatability period. 
            Nd : int
                Nodal days per ground-track repeatability period. 
            observations : list[BaseObservation]
                List of single observation objects to be combined
            )doc")
        .def("get_observations", &MultiObservation::get_observations, R"doc(
            
            Getter for the individual observations

            Return
            --------
            list[BaseObservation]
                List of single observation objects that form the combined observation

            )doc");

    auto gps =
        py::class_<GPS, MultiObservation, std::shared_ptr<GPS>>(m, "GPS");

    gps.def(py::init<int, int, int, double, double, double>(), py::arg("l_max"),
            py::arg("Nr"), py::arg("Nd"), py::arg("I"), py::arg("we_0"),
            py::arg("wo_0"));

    auto constellation =
        py::class_<Constellation, MultiObservation,
                   std::shared_ptr<Constellation>>(m, "Constellation", R"doc(

        This class wraps the :attr:`~gfeatpy.observation.MultiObservation` class to handle
        constellations in a simplified way.
        
        .. note::
            This class assumes that the observation type is :attr:`~gfeatpy.observation.Collinear`.
            Future versions will extend the class to any observation type.
        
        )doc");

    constellation
        .def(py::init<int, int, int, std::vector<double>, double,
                      LongitudePolicy>(),
             py::arg("l_max"), py::arg("Nr"), py::arg("Nd"), py::arg("I"),
             py::arg("rho_0"),
             py::arg("longitude_policy") = LongitudePolicy::INTERLEAVING,
             R"doc(
            
            Parameters
            -----------
            l_max : int
                Cut-off degree.
            Nr : int
                Number of orbital revolutions per ground-track repeatability period. 
            Nd : int
                Nodal days per ground-track repeatability period. 
            I : list[float]
                List of inclinations of the different repeating orbits [rad].
            rho_0 : float
                Nominal intersatellite distance [m].
            longitude_policy : LongitudePolicy
                Longitude separation policy.
             
            )doc")
        .def("set_observation_error", &Constellation::set_observation_error,
             py::arg("range_instrument_asd"), py::arg("accelerometer_asd"),
             R"doc(
            
            An important error contribution to the range observation error is also
            the error in the accelerometer. This error propagates to the range observation.
            This function accounts for this propagation internally.

            
            Parameters
            -----------
            range_instrument_asd : Callable[[float], float]
                Ranging instrument amplitude spectral density for an input frequency.
            accelerometer_asd : Callable[[float], float]
                Accelerometer amplitude spectral density for an input frequency.

            )doc");
}
