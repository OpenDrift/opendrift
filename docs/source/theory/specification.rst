Specifications and background requirements
===========================================

**The goal is to develop a generic, open source trajectory model for any kind of particles (Lagrangian Elements) drifting in the ocean (or even in atmosphere).**


Some desired features are:
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
  * Generic module for **seeding of particles** (point source, line, cone, shapefile, satellite image initialisation etc)
* **Robust**
  * using priority list of driver data in case a source is missing, or if the particles exits the coverage in space or time. Using climatology or constant values as last option.
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
  * platform independent and easy (no) installation
    * based on as few (standard) packages as possible:
      * scipy for scientific functions.
      * matplotlib for plotting functions.
      * basemap for GSHHS coastline for plots and for checking if particles should be stranded.
  * clean code, following the [PEP8](http://www.python.org/dev/peps/pep-0008/) and [PEP257](https://www.python.org/dev/peps/pep-0257/) style guides.

***

See presentations from NOAA regarding GNOME/ADIOS models, with several similar ideas:

https://ceprofs.civil.tamu.edu/ssocolofsky/nfm/Downloads/Workshop_2012/Plenary_Gnome.ppt.pdf
http://www.youtube.com/watch?v=3pJlgQRn1jY
