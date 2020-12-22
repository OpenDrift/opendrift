Using OpenDrift in a container
=====================================

In this tutorial we will walk through building a Docker container with opendrift,
and then interacting with it. If desired, you can use this same container to
create a [Singularity image](https://www.sylabs.io/guides/3.0/user-guide/) that is usable on
a shared resource cluster. We will do the latter by pulling the image on Docker Hub
directly into a Singularity container. This folder contains two Dockerfiles:

 - `Dockerfile <https://github.com/OpenDrift/opendrift/blob/master/Dockerfile>`_: builds a vanilla container with opendrift, Python 3

Pull
----

These images are available on Docker Hub as `opendrift/opendrift <https://hub.docker.com/r/opendrift/opendrift>`_, with instructions provided
below for building locally. You should look at the "tags" tab of each to determine the
version of Open Drift and Python that you are interested in. For example:

 - **opendrift/opendrift:latest** refers to OpenDrift (latest version) with Python 3
 - **opendrift/opendrift:v1.0.7-** refers to OpenDrift (version 1.0.7) with Python 3

If you want to pull a particular version:

.. code::

   $ docker pull opendrift/opendrift:1.0.7

Build
------

If needed, you can develop locally! You should first `install Docker <https://docs.docker.com/install/>`_ so that you
can build images on your host. To build the image, after cloning the repository,
from the base of the repo issue this command, where `opendrift/opendrift` refers
to the name of the container (on Docker Hub it coincides with a `<username><reponame>`
You should select the name of the file that you want to build (e.g., docker/Dockerfile)
and then build as follows:

.. code::

   $ docker build -f docker/Dockerfile -t opendrift/opendrift .

Usage
------

We will show usage examples for the base container, and you can change the name if you
desire to use a different container or version (the tag). Likely you will want to interact with the software, and you can do this via
an interactive (i) terminal (t) session:

.. code::

   $ docker run -it opendrift/opendrift


Singularity
-----------

We will now walk through dumping the Docker layers (you might know Docker images are composed
of layers, like a cake!) into a Singularity container. Why would you want a Singularity container? It doesn't have the security issues of Docker, and could be used on a shared resource.


Installation
++++++++++++
You should first `install Singularity <https://singularityware.github.io/install-linux>`_ so that you can build images on your host. If you use a Mac, you will need to install Singularity in a virtual machine like Vagrant. Singularity is going to allow us to interact exactly the same, but with an image that we can use on our shared resource. The biggest difference is that a Singularity image is a read only single file (a format called `squashfs <https://en.wikipedia.org/wiki/SquashFS>`_ so it is compressed) that we can physically move around and execute like a script.
Unlike Docker images that are assembled from layers and the whole thing is sort of mysterious to the
average user, your Singularity container is a single file that can sit on your desktop.
You can read more about Singularity `here <https://singularityware.github.io>`_.

Pull the Container
++++++++++++++++++
Wherever you are working, the image layers (used to create the Singularity container) will be pulled by default
to the Singularity default cache, which is `$HOME/.singularity`. If there is absolutely any chance that your
`$HOME` has a quota (e.g., on a shared resource) then you should define the environment variable `SINGULARITY_CACHEDIR`
to be somewhere that you **do** have room. For example, it may be the folder defined at the environment variable `$SCRATCH` on your shared resource. You might do something like this before using Singularity:

.. code::

   $ export SINGULARITY_CACHEDIR=$SCRATCH/.singularity

Next, we will pull the Docker Image for OpenDrift from `Docker Hub <https://hub.docker.com/opendrift/opendrift/>`_
directly into a Singularity container. You actually don't need to be a superuser (root / sudo) to do this.

.. code::

   $ singularity pull --name opendrift.sif docker://opendrift/opendrift

Using the image
+++++++++++++++
Let's shell inside the image to interact with the software! I'm not sure how OpenDrift
is used, but I can show you where it is. First, shell inside to explore!


.. code::

   $ singularity shell opendrift.sif


You can also just run the image to get an interactive shell, with working directory `/code`
where the repository was cloned. Here is an example of running directly (without build) to
pull the latest tag:

.. code::

   $ docker run -it opendrift/opendrift

   Unable to find image 'opendrift/opendrift:latest' locally
   latest: Pulling from opendrift/opendrift
   cc1a78bfd46b: Already exists
   bad124d5894e: Pull complete
   ab2b0b173074: Pull complete
   018d53043894: Pull complete
   4987762b1e47: Pull complete
   d58a7f3e3615: Pull complete
   86f53a067a28: Pull complete
   4c17ec80ca72: Pull complete
   aae597ea9e38: Pull complete
   Digest: sha256:33807a79ced6ca9c0960bd942e9d12381c7f1066feb75c5c6992ae5b8802f94c
   Status: Downloaded newer image for opendrift/opendrift:latest
   (base) root@d7ca5fe730b8:/code# python
   Python 2.7.15 |Anaconda, Inc.| (default, May  1 2018, 23:32:55)
   [GCC 7.2.0] on linux2
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import opendrift
   >>> opendrift.__version__
   '1.0.4'

To execute a command to the container from the outside (on the host without shelling
inside) you can use exec:


.. code::

   [vsochat@sh-08-37 ~]$ singularity exec opendrift.sif python myscript.py

The opendrift software (this repository) can be found at `/code/opendrift` in the container.
Note that the creators used / more robustly tested the Singularity container on the
`Sherlock cluster <https://www.sherlock.stanford.edu/docs/>`_.
