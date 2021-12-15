Theory
======

.. toctree::
   :titlesonly:

   lagrangian
   eulerian
   specification
   data_model

Drift in the Ocean
==================

Methods of simulation
---------------------

There are two primary methods of simulating drift in the ocean (as well as in many other fields): *Lagrangian* and *Eulerian*. Simply stated: Lagrangian simulation uses the drifting particle as reference, while Eulerian simulation uses the volume that drift happens through as reference.

In OpenDrift the Lagrangian method is used to simulate particles. This methods requires that you to use a sufficiently large number of particles in order to get a meaningful simulation result.
