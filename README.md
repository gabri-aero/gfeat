# GFEAT


[![PyPI version](https://img.shields.io/pypi/v/gfeatpy)](https://pypi.org/project/gfeatpy/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen)](https://gabri-aero.github.io/gfeat/)
![Linux Only](https://img.shields.io/badge/Linux-Only-red?logo=linux&logoColor=black)






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

The library is available in the Python Package Index. 
```{bash}
pip install gfeatpy
```

### Future Improvements

- Add and improve unit tests for all components.
- Extend Constellation class to any observation.
- Include optimization examples with PyGMO.
- Extend to eccentric orbits.

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

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Balmino1996" class="csl-entry">

Balmino, G., E. Schrama, and N. Sneeuw. 1996.
“<span class="nocase">Compatibility of first-order circular orbit
perturbations theories; consequences for cross-track inclination
functions</span>.” *Journal of Geodesy* 70 (9): 554–61.
<https://doi.org/10.1007/bf00867863>.

</div>

<div id="ref-Jekeli1981" class="csl-entry">

Jekeli, Christopher. 1981. “<span class="nocase">Alternative Methods to
Smooth the Earth’s Gravity Field</span>.” Report 327. Columbus, Ohio:
Department of Geodetic Science; Surveying, The Ohio State University.

</div>

<div id="ref-Kaula1966" class="csl-entry">

Kaula, William M. 1966. *<span class="nocase">Theory of Satellite
Geodesy: Applications of Satellites to Geodesy</span>*. Blaisdell
Publishing Company.

</div>

<div id="ref-Schrama1990" class="csl-entry">

Schrama, E. J. O. 1990. “<span class="nocase">Gravity field error
analysis: Applications of GPS receivers and gradiometers on low orbiting
platforms</span>.” NASA.

</div>

<div id="ref-Schrama1991" class="csl-entry">

———. 1991a. “<span class="nocase">Error Propagation and Correlation
Analysis of Covariance Matrices</span>.” Goddard Space Flight Center;
Personal communication.

</div>

<div id="ref-Schrama1991a" class="csl-entry">

———. 1991b. “<span class="nocase">Gravity field error analysis:
Applications of global positioning system receivers and gradiometers on
low orbiting platforms</span>.” *Journal of Geophysical Research: Solid
Earth* 96 (B12): 20041–51. <https://doi.org/10.1029/91jb01972>.

</div>

<div id="ref-Sneeuw1994" class="csl-entry">

Sneeuw, Nico. 1994. “<span class="nocase">Global spherical harmonic
analysis by least-squares and numerical quadrature methods in historical
perspective</span>.” *Geophysical Journal International* 118 (3):
707–16. <https://doi.org/10.1111/j.1365-246x.1994.tb03995.x>.

</div>

<div id="ref-Sneeuw2000" class="csl-entry">

Sneeuw, Nicolaas. 2000. “<span class="nocase">A Semi-Analytical Approach
to Gravity Field Analysis from Satellite Observations</span>.” PhD
thesis, Institut für Astronomische und Physikalische Geodäsie, TU
München.

</div>

<div id="ref-Visser2012" class="csl-entry">

Visser, P. N. A. M., E. J. O. Schrama, N. Sneeuw, and M. Weigelt. 2012.
“<span class="nocase">Dependency of Resolvable Gravitational Spatial
Resolution on Space-Borne Observation Techniques</span>.” In
*<span class="nocase">Geodesy for Planet Earth</span>*, edited by Steve
Kenyon, Maria Christina Pacino, and Urs Marti, 373–79. Berlin,
Heidelberg: Springer Berlin Heidelberg.

</div>

<div id="ref-Wahr1998" class="csl-entry">

Wahr, John, Mery Molenaar, and Frank Bryan. 1998.
“<span class="nocase">Time variability of the Earth’s gravity field:
Hydrological and oceanic effects and their possible detection using
GRACE</span>.” *Journal of Geophysical Research: Solid Earth* 103 (B12):
30205–29. <https://doi.org/10.1029/98jb02844>.

</div>

</div>



