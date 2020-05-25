Installing OpenDrift
=============================================

OpenDrift is usually installed as a source package or as a conda package.

.. note::
   If you expect to make changes to the the OpenDrift code (e.g. models) it is
   more convenient to install from :ref:`source <source_install>`.

.. warning::

   OpenDrift has now transitioned to *Python 3*. New versions of OpenDrift are unlikely to work with Python 2.
   Python 2 reached end of line in January, 2020. If you need to use Python 2 try using an older version of OpenDrift.


.. _miniconda_install:

Miniconda (recommended)
-----------------------

1. Install `miniconda3 <https://docs.conda.io/en/latest/miniconda.html>`_
2. Set up a *Python 3* environment for opendrift

.. code-block:: bash

   $ conda config --add channels conda-forge  # recommended, but not necessary
   $ conda create -n opendrift python=3
   $ conda activate opendrift

.. note::
   Windows users might have to specify `python=3.7.4`.

3. Install OpenDrift and dependencies

.. code-block:: bash

  $ conda install -c opendrift -c conda-forge -c noaa-orr-erd opendrift

.. _source_install:

From source
-----------

Using Miniconda (recommended)
+++++++++++++++++++++++++++++

.. note::
   There are two main ways of installing OpenDrift from source in a conda environment. If you already have OpenDrift installed and wish to
   use the development version, the first is preferred. If you wish to set up to run from source, the second method is probably more suitable (and
   corresponds to how OpenDrift is setup when developed).

Alternative 1: Already installed
_________________

1. Follow the instructions in :ref:`miniconda_install`.
2. Clone OpenDrift:

.. code-block:: bash

   $ git clone https://github.com/OpenDrift/opendrift.git
   $ cd opendrift/

3. Install as editable:

.. code-block:: bash

   $ pip install -e .

Alternative 2: Fresh development environment
_____________________________

1. Install `miniconda3 <https://docs.conda.io/en/latest/miniconda.html>`_
2. Clone OpenDrift:

.. code-block:: bash

   $ git clone https://github.com/OpenDrift/opendrift.git
   $ cd opendrift/

3. Create environment with the required dependencies and install OpenDrift

.. code-block:: bash

  $ conda config --add channels conda-forge  # recommended, but not necessary
  $ conda env create -f environment.yml
  $ conda activate opendrift
  $ pip install -e .

This installs the OpenDrift package as an editable package. You can therefore directly make changes to the repository or fetch the newest changes with :code:`git pull`. You do not need to add OpenDrift to PYTHONPATH as long as you have the :code:`opendrift` environment activated.


