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

OpenDrift is usually installed as a source package or as a conda package. It is
also possible to install via `pypi.org <http://pypi.org>`_.

Miniconda
-----------

1. Install `miniconda3 <https://docs.conda.io/en/latest/miniconda.html>`_
2. Set up a *Python 3* environment for opendrift

.. code-block:: bash

   $ conda env create -n opendrift python=3
   $ conda activate opendrift

3. Install opendrift

.. code-block:: bash

  $ conda install -c opendrift

.. If you exepect to make changes to the the OpenDrift code (e.g. models) it is
.. probably more convenient to install from source: :ref:`install_source`.
.. See :doc:`install/pip` or :doc:`install/source`.


Contents
===============

.. toctree::
   :maxdepth: 3
   :glob:

   install
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
