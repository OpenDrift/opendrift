# OpenDrift in a Container 

In this tutorial we will walk through building a Docker container with opendrift,
and then interacting with it. If desired, you can use this same container to 
create a [Singularity image](https://singularityware.github.io) that is usable on
a shared resource cluster. We will do the latter by pulling the image on Docker Hub
directly into a Singularity container. This folder contains two Dockerfiles:

 - [Dockerfile](Dockerfile): builds a vanilla container with just opendrift
 - [Dockerfile.oil](Dockerfile.oil) builds the same container with opendrift + NOAA oillibrary for oil simulations.

## Docker 

You should first [install Docker](https://docs.docker.com/install/) so that you 
can build images on your host. 

### Build
To build the image, after cloning the repository,
from the base of the repo issue this command, where `vanessa/opendrift` refers
to the name of the container (on Docker Hub it coincides with a `<username><reponame>`

```bash
## This is how the container was build
$ docker build -t vanessa/opendrift .

## NOAA OILLibrary
$ docker build -f Dockerfile.oil -t vanessa/opendrift-oil .
```

### Pre-built
If you don't want to take this step, the containers are provided on Dockerhub, ([opendrift](https://hub.docker.com/r/vanessa/opendrift/),[opendrift-oil](https://hub.docker.com/r/vanessa/opendrift-oil/)) and you can just run or pull either of them. 
These containers were built locally and pushed directly,
and currently there are the following versions:

 - **vanessa/opendrift:1.0.3** refers to version 1.0.3
 - **vanessa/opendrift:1.0.4** refers to version 1.0.4
 - **vanessa/opendrift:latest** also refers to version 1.0.4
 - **vanessa/opendrift-oil:1.0.4** refers to version 1.0.4
 - **vanessa/opendrift-oil:latest** also refers to version 1.0.4

This means that you can skip the build step above! If you want to pull a particular version:

```bash
$ docker pull vanessa/opendrift:1.0.4
$ docker pull vanessa/opendrift-oil:1.0.4
```

And if there is a newer/older version you'd like to request built, if the OpenDrift maintainers
are not providing it, please reach out to [@vsoch](https://www.github.com/vsoch) and she will be happy to build and push
an updated version for you.

### Usage
We will show usage examples for the base container, and you can change the name if you
desire to use a different container or version (the tag). Likely you will want to interact with the software, and you can do this via 
an interactive (i) terminal (t) session:

```bash
$ docker run -it vanessa/opendrift
```

The original container recipe (Dockerfile) and documentation are provided [here](https://github.com/researchapps/sherlock/tree/master/opendrift) if you need to ask for help from the creator.


## Singularity

We will now walk through dumping the Docker layers (you might know Docker images are composed
of layers, like a cake!) into a Singularity container. Why would you want a Singularity container? It doesn't have the security issues of Docker, and could be used on a shared resource.


## Installation
You should first [install Singularity](https://singularityware.github.io/install-linux) so that you can build images on your host. If you use a Mac, you will need to install Singularity in a virtual machine like Vagrant. Singularity is going to allow us to interact exactly the same, but with an image that we can use on our shared resource. The biggest difference is that a Singularity image is a read only single file (a format called [squashfs](https://en.wikipedia.org/wiki/SquashFS) so it is compressed) that we can physically move around and execute like a script.
Unlike Docker images that are assembled from layers and the whole thing is sort of mysterious to the 
average user, your Singularity container is a single file that can sit on your desktop. 
You can read more about Singularity [here](https://singularityware.github.io). 

### Pull the Container
Wherever you are working, the image layers (used to create the Singularity container) will be pulled by default
to the Singularity default cache, which is `$HOME/.singularity`. If there is absolutely any chance that your
`$HOME` has a quota (e.g., on a shared resource) then you should define the environment variable `SINGULARITY_CACHEDIR`
to be somewhere that you **do** have room. For example, it may be the folder defined at the environment variable `$SCRATCH` on your shared resource. You might do something like this before using Singularity:

```bash
$ export SINGULARITY_CACHEDIR=$SCRATCH/.singularity
```

Next, we will pull the Docker Image for OpenDrift from [Docker Hub](https://hub.docker.com/vanessa/opendrift/)
directly into a Singularity container. You actually don't need to be a superuser (root / sudo) to do this.

```bash
$ singularity pull --name opendrift.simg docker://vanessa/opendrift
```

### Using the image
Let's shell inside the image to interact with the software! I'm not sure how OpenDrift
is used, but I can show you where it is. First, shell inside to explore!


```bash
$ singularity shell opendrift.simg
python
Python 2.7.15 |Anaconda, Inc.| (default, May  1 2018, 23:32:55) 
[GCC 7.2.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import opendrift
>>> 
```


You can also just run the image to get an interactive shell, with working directory `/code`
where the repository was cloned. Here is an example of running directly (without build) to
pull the latest tag:

```bash
$ docker run -it vanessa/opendrift
Unable to find image 'vanessa/opendrift:latest' locally
latest: Pulling from vanessa/opendrift
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
Status: Downloaded newer image for vanessa/opendrift:latest
(base) root@d7ca5fe730b8:/code# python
Python 2.7.15 |Anaconda, Inc.| (default, May  1 2018, 23:32:55) 
[GCC 7.2.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import opendrift
>>> opendrift.__version__
'1.0.4'
```

To execute a command to the container from the outside (on the host without shelling
inside) you can use exec:


```bash
[vsochat@sh-08-37 ~]$ singularity exec opendrift.simg python myscript.py
```

The opendrift software (this repository) can be found at `/code/opendrift` in the container.
Note that the creators used / more robustly tested the Singularity container on the 
[Sherlock cluster](https://www.sherlock.stanford.edu/docs/).
