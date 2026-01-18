.. GFEAT documentation master file, created by
   sphinx-quickstart on Mon Jul  7 18:19:27 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Examples
========

Gravity field data processing examples
--------------------------------------

These examples show how to load and process different types of gravity field data.

.. nbgallery::
   
   notebooks/gravity-field-data/static.ipynb
   notebooks/gravity-field-data/monthly.ipynb
   notebooks/gravity-field-data/normals.ipynb
   notebooks/gravity-field-data/aod1b.ipynb


Gravity field error analysis
-------------------------------------

Here, the utilities of the observation module are showcased. It consists of an analytical gravity field error methodology 
based on the lumped coefficients' theory (Colombo, 1984).

Single satellite pair
~~~~~~~~~~~~~~~~~~~~~~

These examples investigate the gravity field recovery capabilities of inter-satellite range observations for a single-pair 
collinear formation. Some fundamental properties of the analytical methodology are also described here.

.. nbgallery::
   
   notebooks/single-pair/sensitivity.ipynb
   notebooks/single-pair/matrix-structure.ipynb
   notebooks/single-pair/regularized.ipynb
   notebooks/single-pair/grace.ipynb
   notebooks/single-pair/grace-resonant.ipynb
   notebooks/single-pair/ground-track-analysis.ipynb

Multi-pair constellations
~~~~~~~~~~~~~~~~~~~~~~~~~

These examples analyse the capabilities of multi-pair satellite constellations. Only inter-satellite ranging in a collinear
formation is considered. Crucial mission design aspects are observed, such as the impact of non-polar inclinations, ground-track 
policies or the inherent anisotropy of the gravity field solution for such formations.

.. nbgallery::
   
   notebooks/multiple-pair/bender.ipynb
   notebooks/multiple-pair/longitude-policy.ipynb

