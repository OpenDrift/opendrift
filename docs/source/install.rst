Installing OpenDrift
=============================================

OpenDrift is usually installed as a source package or as a conda package. It is
also possible to install via `pypi.org <http://pypi.org>`_.

.. warning::

   OpenDrift has now transitioned to *Python 3*. New versions of OpenDrift are unlikely to work with Python 2.
   Python 2 reaches end of line in January, 2020. If you need to use Python 2 try using an older version of OpenDrift.


Miniconda
-----------

1. Install `miniconda3 <https://docs.conda.io/en/latest/miniconda.html>`_
2. Set up a *Python 3* environment for opendrift

.. code-block:: bash

   $ conda create -n opendrift python=3
   $ conda activate opendrift

3. Install OpenDrift and dependencies

.. code-block:: bash

  $ conda install -c opendrift -c conda-forge -c noaa-orr-erd opendrift


.. warning::
   If you exepect to make changes to the the OpenDrift code (e.g. models) it is
   probably more convenient to install from source.


Using pip
----------

1. Consider creating an exclusive virtual environment for OpenDrift.
2. Install dependencies for the OilLibrary:

.. code-block:: bash

  $ pip install gitpython sqlalchemy zope.sqlalchemy==1.1 transaction numpy scipy https://github.com/NOAA-ORR-ERD/PyNUCOS/archive/v2.5.5.tar.gz awesome-slugify

3. Install the OpenDrift version of the OilLibrary

.. code-block:: bash

  $ pip install git+https://github.com/OpenDrift/OilLibrary

4. Install OpenDrift

.. code-block:: bash

  $ pip install opendrift

From source
-----------

1. Consider creating an exclusive virtual environment for OpenDrift.

2. Check out source

.. code-block:: bash

  $ git clone https://github.com/OpenDrift/opendrift.git

3. Install dependencies for the OilLibrary:

.. code-block:: bash

  $ pip install gitpython sqlalchemy zope.sqlalchemy==1.1 transaction numpy scipy https://github.com/NOAA-ORR-ERD/PyNUCOS/archive/v2.5.5.tar.gz awesome-slugify

4. Install the OpenDrift version of the OilLibrary

.. code-block:: bash

  $ pip install git+https://github.com/OpenDrift/OilLibrary

5. Install OpenDrift dependencies

.. code-block:: bash

  $ pip install -r requirements.txt

6. Install OpenDrift as an editable package

.. code-block:: bash

  $ pip install -e .

