# GFEAT


[![PyPI version](https://img.shields.io/pypi/v/gfeatpy)](https://pypi.org/project/gfeatpy/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen)](https://gabri-aero.github.io/gfeat/)


> <b>G</b>ravity <b>F</b>ield <b>E</b>rror <b>A</b>nalysis <b>T</b>ool (GFEAT) is a C++ library that employs an analytical method to estimate gravity field error from different orbital configurations and observations.

### Features

- Estimation of gravity field error.
    - Error propagation from observations PSD to the Stokes' coefficients.
    - Several observation types supported: gravity potential, radial, along-track and cross-track displacements, inter-satellite range...
    - Combination of multiple observation types: 3D GPS position and multi-pair constellations.
- Processing of gravity field data.
    - Gravity field I/O handling.
    - Error covariance handling from normal equations.
    - Synthesis with different functionals: EWH, gravity anomalies, geoid heights.
    - Covariance synthesis also supported.
    - FFT employed for efficient computation.


### Installation

The library is available in the Python Package Index. It is recommended to create a conda environment first.
```{bash}
conda create env --name gfeat_env
```
Then, install the library.
```{bash}
pip install gfeatpy
```

### Future Improvements

- Add and improve unit tests for all components.
- Extend Constellation class to any observation.
- Include range-rate, range-acceleration.
- Include optimization examples with PyGMO.

### C++ development

To use the library natively in C++ is currently not recommended given the limitations for output data as well as the easier to use plotting tools available in Python. However, a good starting point are the C++ unit tests, under the `test` folder.

Some very descriptive examples can be found in [TestConstellation.cpp](https://github.com/gabri-aero/gfeat/blob/main/tests/observation/TestConstellation.cpp) and [TestCollinear.cpp](https://github.com/gabri-aero/gfeat/blob/main/tests/observation/TestCollinear.cpp). You can build and run all the tests like this.
```{make}
mkdir build
cd build
cmake ..
make
ctest
```


### References




