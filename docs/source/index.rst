.. toctree::
   :hidden:

   self

OpenDrift
=====================================

OpenDrift is a software for modeling the trajectories and fate of objects or
substances drifting in the ocean, or even in the atmosphere.

OpenDrift is open source, and is programmed in Python. As the software is very
generic, it is rather a "framework" than a "trajectory model" in the
traditional sense. Trajectory models for specific purposes (e.g. oil drift,
search and rescue, larvae drift etc) may reuse all common functionality from
the core model, and need only implement a Python Class describing the
purpose-specific processes (physics/biology etc). See :doc:`theory/specification`
and :doc:`theory/data_model` for more detailed information.



Installation
=============


Contents
===============

.. toctree::
   :maxdepth: 3
   :glob:

   theory/specification
   theory/data_model
   models
   autoapi/index


Indices and tables
==================

* :ref:`genindex`
.. * :ref:`modindex`
* :ref:`search`
* :doc:`autoapi/index`
