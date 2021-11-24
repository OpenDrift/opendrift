Specifications and background requirements
===========================================

**OpenDrift is designed to be a generic, open source trajectory model for any kind of particles (Lagrangian Elements) drifting in the ocean (or even in atmosphere).**


Some design criteria are:

* **Modular**

  * Can read input (driver) data (currents, wind, waves etc) from any source/model in any projection and file format (including OPeNDAP/Thredds)

    * reads data only for the relevant locations/times where the active drifting objects are located
    * should be easy to add “readers” for new sources/models
    * parameters referred to by their CF standard-names, for consistency

  * A **core module** will contain everything common to all types of particles

    * fetch driver data (wind, wave, current, water temp, air temp, diffusivity, etc)
    * check beaching/stranding
    * export output, with possibility to import for continued simulation
    * switch between "trajectory view" (history for each particle) and "field view" (particle positions/properties at given time)
    * plot figure (trajectories, or spatial distribution at given time; density maps, convex hull, ...)

  * **Separate modules** with functions to advect and update the properties of specific types of particles:

    * oil (based on OD3D)
    * rigid objects (based on Leeway)
    * plankton, fish eggs, larvae (based on Ladim)
  * Generic methods for **seeding of particles** (point source, line, cone, shapefile, satellite image initialisation etc)
 
* **Robust**

  * using priority list of driver data in case a source is missing, or if the particles leaves the domain in space or time. Using climatology or constant values as last option.

  * gracious error handling (try-except)
 
* **Fast** (enough)

  * Using vectorized code (performance similar to c and fortran)
 
* **Flexible**

  * Can be run for at any location worldwide where driver data are available
  * Option to run both forwards and backwards in time
  * Modules and functions/parameterisations can be modified easily by the user, e.g. for sensitivity tests and scientific studies
 
* **Stochastic**

  * Take into account uncertainty:

    * ensemble input
    * multi-driver ensemble
    * error, standard deviation

  * Yet reproducible, using pseudo-random numbers, allowing sensitivity tests (various parameterisations, driver data sources etc)
 
* **Programmed in Python**

  * with optional possibility to implement specific modules og bottlenecks in other laguages, e.g. c or fortran
  * platform independent and simple installation

    * based on as few external packages as possible

  * clean code, following the `PEP8 <https://www.python.org/dev/peps/pep-0008>`_ and `PEP257 <https://www.python.org/dev/peps/pep-0257/>`_ style guides.
