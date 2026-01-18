Planet object
===============

.. automodule:: gfeatpy
   :members:
   :show-inheritance:

.. code-block:: python

   from gfeatpy import planet    # mind the lower case 

   # Asign central planetary body constants to Mars
   planet.mu = 42.83e12
   planet.ae = 3396.2e3
   planet.C20 = -1960.45e-6
   planet.theta_dot = 7.088e-5
   planet.rho_e = 2582



Gravity module
===============

The gravity module provides different tools to handle gravity field data.


.. automodule:: gfeatpy.gravity
   :no-members:

.. autoclass:: gfeatpy.gravity.AOD1B
   :members:

.. py:enum:: AOD1BType

   Enumeration that defines the Stokes coefficients set to be retrieved from AOD1B data.

   .. py:attribute:: ATM
   
      Difference between vertically integrated density of the atmosphere and the corresponding mean field.

   .. py:attribute:: OCN
   
      Difference between the water column contribution to ocean bottom pressure and the corresponding mean field.

   .. py:attribute:: GLO
   
      Sum of ATM and OCN mass anomalies.

   .. py:attribute:: OBA
   
      Sum of the water column contribution to the ocean bottom pressure anomalies (OCN) 
      and the atmospheric contribution to ocean bottom pressure anomalies.

.. autoclass:: gfeatpy.gravity.BaseFunctional

   .. inheritance-diagram::
      gfeatpy.gravity.GeoidHeight
      gfeatpy.gravity.GravityAnomaly
      gfeatpy.gravity.EquivalentWaterHeight
      :top-classes: gfeatpy.gravity.BaseFunctional
      :parts: 1

.. autoclass:: gfeatpy.gravity.DateTime
   :members:

.. autoclass:: gfeatpy.gravity.EquivalentWaterHeight
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.gravity.GeoidHeight
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.gravity.GravityAnomaly
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.gravity.GravityField
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.gravity.SphericalHarmonics
   :members:

.. autoclass:: gfeatpy.gravity.SphericalHarmonicsCovariance
   :members:

Observation module
==================

The observation module handles the computation of the spectral observations. For this purpose, different functionalities
are defined in abstract objects and inherited towards the actual usable observation classes.

.. inheritance-diagram:: 
   gfeatpy.observation.AbstractKiteSystem
   gfeatpy.observation.BaseObservation
   gfeatpy.observation.MultiObservation
   gfeatpy.observation.Radial
   gfeatpy.observation.CrossTrack
   gfeatpy.observation.AlongTrack
   gfeatpy.observation.Potential
   gfeatpy.observation.Range
   gfeatpy.observation.Constellation
   gfeatpy.observation.GPS
   :top-classes: gfeatpy.observation.AbstractKiteSystem
   :parts: 1

The observations require the definition of different orbital elements which are relative to the Earth frame as depicted below.

.. image:: _static/orbital-elements.jpg
   :alt: Orbital elements figure
   :width: 400px

.. automodule:: gfeatpy.observation
   :no-members:

.. autoclass:: gfeatpy.observation.AbstractKiteSystem
   :members:

.. autoclass:: gfeatpy.observation.AlongTrack
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.observation.BaseObservation
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.observation.Constellation
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.observation.GPS
   :members:
   :show-inheritance:

.. py:enum:: LongitudePolicy

   Enumeration that defines the two built-in longitude separation policies for the :attr:`~gfeatpy.observation.Constellation`
   class.

   .. py:attribute:: INTERLEAVING
   
      The ground-tracks are perfectly interleaved at the equator.

   .. py:attribute:: OVERLAPPING
   
      The ground-tracks overlap at the equator.

.. autoclass:: gfeatpy.observation.MultiObservation
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.observation.Potential
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.observation.Radial
   :members:
   :show-inheritance:

.. autoclass:: gfeatpy.observation.Range
   :members:
   :show-inheritance:

Plotting module
===============

.. automodule:: gfeatpy.plotting
   :members:
   :show-inheritance:

Utils module
===============

.. automodule:: gfeatpy.utils
   :members:
   :show-inheritance:
