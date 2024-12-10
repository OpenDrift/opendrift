Installing OpenDrift
=============================================

Alternative 1: Using Mambaforge (or Miniforge) and Git (recommended)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

1. Install `miniforge <https://conda-forge.org/miniforge/>`_
2. Clone OpenDrift:

.. code-block:: bash

   $ git clone https://github.com/OpenDrift/opendrift.git
   $ cd opendrift/

3. Create environment with the required dependencies and install OpenDrift

.. code-block:: bash

  $ mamba env create -f environment.yml
  $ conda activate opendrift
  $ pip install --no-deps -e .

This installs the OpenDrift package as an editable package. You can therefore directly make changes to the repository or fetch the newest changes with :code:`git pull`. You do not need to add OpenDrift to PYTHONPATH as long as you have the :code:`opendrift` environment activated.

Occasionally, new dependencies will be added to environment.yml. Then the local conda environment can be updated with:

.. code-block:: bash

  $ mamba env update -f environment.yml

Alternative 2: Using Mambaforge (or Miniforge)
++++++++++++++++++++++++++++++++++++++++++++++

1. Install `miniforge <https://conda-forge.org/miniforge/>`_
2. Set up a *Python 3* environment for opendrift and install opendrift

.. code-block:: bash

   $ mamba create -n opendrift python=3 opendrift
   $ conda activate opendrift

.. _source_install:

Modify install to use Git version
---------------------------------

If you later want to edit the OpenDrift source code, or be able to update from repository with `git pull`, the following two steps are necessary. This yields the same as Alternative 1.

3. Clone OpenDrift:

.. code-block:: bash

   $ git clone https://github.com/OpenDrift/opendrift.git
   $ cd opendrift/

4. Remove original package:

.. code-block:: bash

   $ mamba remove --force opendrift

5. Install as editable:

.. code-block:: bash

   $ pip install --no-deps -e .

Building and using the Docker image
+++++++++++++++++++++++++++++++++++

It is also possible to run OpenDrift in a `docker <https://www.docker.com/>`_
container. You can either use the `pre-built images (each release or latest
master) <https://hub.docker.com/r/opendrift/opendrift>`_:

.. code-block:: bash

  $ docker run -it --rm opendrift/opendrift:latest

or if you want to re-build the image using your local version of opendrift, check out opendrift and make your changes, then:

.. code-block:: bash

   $ git clone https://github.com/OpenDrift/opendrift.git
   $ cd opendrift/

re-build the image:

.. code-block:: bash

  $ docker build -t opendrift .

and run it:

.. code-block:: bash

  $ docker run -it --rm opendrift

