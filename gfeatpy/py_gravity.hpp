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
            different gravity field functionals. For this purpose, it contains a function
            ``degree_dependent_factor`` that is overloaded by the child classes. 

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

    base_functional.def("_degree_dependent_factor",
                        &BaseFunctional::degree_dependent_factor,
                        py::arg("degree"));

    base_functional.def_readonly("_common_factor",
                                 &BaseFunctional::common_factor);
    gravity_anomaly.def(py::init<>(),
                        R"doc(
            Constructor for GravityAnomaly.
                
            Example
            --------
            Create a BaseFunctional instance for gravity anomalies::

                gravity_anomaly = GravityAnomaly()
    )doc");
    geoid_height.def(py::init<>(),
                     R"doc(
            Constructor for GeoidHeight.
                
            Example
            --------
            Create a BaseFunctional instance for geoid heights::

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

    auto sh = py::class_<SphericalHarmonics>(m, "SphericalHarmonics", R"doc(
            This class stores Spherical Harmonics coefficients data and provides utilities to 
            process it. Data is loaded through auxiliar objects such as :attr:`~GravityField` or :attr:`~AOD1B`.
        )doc");

    sh.def(py::init<int>(), py::arg("l_max"), R"doc(
            Contructor for SphericalHarmonics object.
            The constructor only defines the internal storing structure.

            Parameters
            -----------
            l_max : int
                Cut-off degree.
        )doc")
        .def("potential", &SphericalHarmonics::potential, py::arg("r_ecef"),
             R"doc(
            Function that evaluates gravity field perturbing potential in the Earth-Centered-Earth-Fixed (ECEF) frame .

            Parameters
            -----------
            r_ecef : numpy.ndarray
                3D position in the ECEF frame.

            Return
            --------
            float
                Gravity field perturbing potential.
        )doc")
        .def("gravity", &SphericalHarmonics::gravity, py::arg("r_ecef"),
             R"doc(
            Function that evaluates the gravity vector in the Earth-Centered-Earth-Fixed frame (ECEF).

            Parameters
            -----------
            r_ecef : numpy.ndarray
                3D position in the ECEF frame.

            Return
            --------
            numpy.ndarray
                3D gravity vector in the ECEF frame.
        )doc")
        .def("synthesis", &SphericalHarmonics::synthesis, R"doc(
            This function performs Global Spherical Harmonics Synthesis according to the 
            input gravity field functional.
            
            .. math::

                y(\lambda, \phi) = \sum_{l=2}^{L} \sum_{m=0}^{l} f(l) \bar{P}_{lm}(\sin{\phi}) (\bar{C}_{lm} \cos{m\lambda} + \bar{S}_{lm} \sin{m\lambda})

            It makes use of the two-step fold proposed by Rizos (1979) as well as the FFT for efficient computation.

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
                A 2D array with functional values on the lon/lat grid.
             )doc")
        .def("degree_variance", &SphericalHarmonics::degree_variance,
             py::arg("use_sigmas"),
             R"doc(
            This function computes degree variance spectrum from the values 
            of the SH coefficients.

            .. math::
                \sigma_l^2 = \sum_{m=0}^{l} \bar{C}_{lm}^2 + \bar{S}_{lm}^2

            Parameters
            -----------
            use_sigmas : bool
                Flag to indicate if the standard deviations are used instead of the coefficients.
            )doc")
        .def("rms_per_coefficient_per_degree",
             &SphericalHarmonics::rms_per_coefficient_per_degree,
             py::arg("use_sigmas"),
             R"doc(
            This function computes RMS per coefficient per degree from the values 
            of the SH coefficients.

            .. math::
                \delta_l = \sqrt{\sum_{m=0}^{l} \frac{\bar{C}_{lm}^2 + \bar{S}_{lm}^2}{2l+1}}

            Parameters
            -----------
            use_sigmas : bool
                Flag to indicate if the standard deviations are used instead of the coefficients.
            )doc")
        .def_readwrite("mu", &SphericalHarmonics::mu, "")
        .def_readwrite("ae", &SphericalHarmonics::ae)
        .def_readwrite("coefficients", &SphericalHarmonics::coefficients,
                       R"doc(
            This property points to the (L+1)x(L+1) dense matrix storing the SH coefficients. The matrix
            structure is as follows.

            .. math::
                \Sigma_{xx} = \begin{bmatrix}
                    \bar{C}_{00} & \bar{S}_{11} & \cdots & \bar{S}_{L1} \\
                    \bar{C}_{10} & \bar{C}_{11} & \cdots & \bar{S}_{L2} \\
                    \vdots & \vdots & \ddots & \vdots \\
                    \bar{C}_{L0} & \bar{C}_{L1} & \cdots & \bar{C}_{LL}
                \end{bmatrix}
            )doc")
        .def_readwrite("sigmas", &SphericalHarmonics::sigmas,
                       R"doc(
            This property points to the (L+1)x(L+1) dense matrix storing the standard deviations of the
            SH coefficients. The matrix structure is equivalent to :py:attr:`~SphericalHarmonics.coefficients`.
            )doc");
    ;
    auto sh_error = py::class_<SphericalHarmonicsCovariance>(
        m, "SphericalHarmonicsCovariance", R"doc(
            This class stores the full covariance matrix from a gravity field solution and provides utilities
            to process it.
        )doc");

    sh_error
        .def(py::init<int>(), py::arg("l_max"),
             R"doc(
            Constructor for SphericalHarmonicsCovariance object. 
            The constructor only defines the inner structure of the matrix. Data can be load through different
            loading functions depending on the file type.

            Parameters
            -----------
            l_max : int
                Cut-off degree.
        )doc")
        .def("from_normal", &SphericalHarmonicsCovariance::from_normal,
             py::arg("filename"), py::arg("scaling_factor"),
             R"doc(
            Function to fill the SH covariance matrix through inversion of the normal matrix as provided
            by an input Sinex file.

            Parameters
            -----------
            filename : string
                Path to Sinex file that contains the normal matrix.
            scaling_factor : float
                Scaling factor to be applied to the normal matrix.         
            )doc")
        .def("synthesis", &SphericalHarmonicsCovariance::synthesis,
             py::arg("n_lon"), py::arg("n_lat"), py::arg("functional"), R"doc(
            This function performs covariance propagation from the SH covariance into the sphere according to the 
            input gravity field functional.
            
            .. math::

                \sigma_{y}(\lambda, \phi) = v P_{xx} v^T

            It makes use of the two-step fold proposed by Rizos (1979) as well as the FFT for efficient computation.

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
        .def("degree_variance", &SphericalHarmonicsCovariance::degree_variance,
             R"doc(
            This function computes degree variance spectrum from the standard deviation 
            of the SH coefficients.

            .. math::
                \sigma_l^2 = \sum_{m=0}^{l} \sigma^2(\bar{C}_{lm}) + \sigma^2(\bar{S}_{lm})
            )doc")
        .def("rms_per_coefficient_per_degree",
             &SphericalHarmonicsCovariance::rms_per_coefficient_per_degree,
             R"doc(
            This function computes RMS per coefficient per degree from the standard deviation 
            of the SH coefficients.

            .. math::
                \delta_l = \sqrt{\sum_{m=0}^{l} \frac{\sigma^2(\bar{C}_{lm}) + \sigma^2(\bar{S}_{lm})}{2l+1}}
            )doc")
        .def_readwrite("Pxx", &SphericalHarmonicsCovariance::Pxx, R"doc(
            Attribute that stores the full covariance matrix.
            )doc");

    auto gravity_field =
        py::class_<GravityField, SphericalHarmonics>(m, "GravityField", R"doc(
            This is a derived class from :attr:`~SphericalHarmonics` that includes functionality
            to load data from .gfc files (see `ICGEM website <https://icgem.gfz-potsdam.de/>`_).
            )doc");

    gravity_field
        .def(py::init<int>(), py::arg("l_max"), R"doc(
            Contructor for GravityField object.

            Parameters
            -----------
            l_max : int
                Cut-off degree.
        )doc")
        .def("load", &GravityField::load, py::arg("filename"), R"doc(
            This function loads gravity field data from a .gfc file into this object.

            Parameters
            -----------
            filename : string
                Path to .gfc file.

            Return
            --------
            GravityField
                Object with loaded gravity field.
        )doc");

    auto datetime = py::class_<DateTime>(m, "DateTime", R"doc(
            Class that defines date and time from calendar format. 
        )doc");

    datetime
        .def(py::init<int, int, int, int, int, int>(), py::arg("year"),
             py::arg("month"), py::arg("day"), py::arg("h") = 0,
             py::arg("m") = 0, py::arg("s") = 0, R"doc(
            
            Constructor for DateTime object.

            Parameters
            -----------
            year : int
            month : int
            day : int
            hours : int
            minutes : int
            seconds : int

             
            )doc")
        .def_readwrite("year", &DateTime::year)
        .def_readwrite("month", &DateTime::month)
        .def_readwrite("day", &DateTime::day)
        .def_readwrite("h", &DateTime::h)
        .def_readwrite("m", &DateTime::m)
        .def_readwrite("s", &DateTime::s);

    auto aod1b_type = py::enum_<AOD1BType>(m, "AOD1BType", R"doc(
            Enumeration that defines the Stokes coefficients set to be retrieved from AOD1B data as defined by the
            AOD1B Product Description Document for Product Release 07
        )doc");

    aod1b_type
        .value("ATM", AOD1BType::ATM, R"doc(
            Difference between vertically integrated density of the atmosphere and the corresponding mean field.
        )doc")
        .value("OCN", AOD1BType::OCN, R"doc(
            Difference between the water column contribution to ocean bottom pressure and the corresponding mean field.
            )doc")
        .value("GLO", AOD1BType::GLO, R"doc(
            Sum of ATM and OCN mass anomalies. 
            )doc")
        .value("OBA", AOD1BType::OBA, R"doc(
            Sum of the water column contribution to the ocean bottom pressure anomalies (OCN) 
            and the atmospheric contribution to ocean bottom pressure anomalies (it is set to zero
            at the continents).
            )doc")
        .export_values();

    auto aod1b = py::class_<AOD1B>(m, "AOD1B");

    aod1b.def(py::init<>())
        .def("load", &AOD1B::load, py::arg("filename"), R"doc(
            Function that loads an AOD1B daily .asc file into the AOD1B class.

            Parameters
            -----------
            filename : string
                Path to AOD1B .asc file.

            Return
            --------
            AOD1B
                AOD1B object with loaded data.
        )doc")
        .def("get", &AOD1B::get, py::arg("datetime"), py::arg("type"), R"doc(
            Function that retrieves the :attr:`~SphericalHarmonics` object associated to the AOD1B
            set indicated by the input parameters.

            Parameters
            -----------
            datetime : DateTime
                Epoch associated to the solution (note that they are provided every 3h).
            type : AOD1BType
                Stokes coefficient set type from the four different possibilities provided.

            Return
            --------
            SphericalHarmonics
                Object containing the associated Stokes coefficients.

            )doc");
}
