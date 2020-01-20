Data model
==========

The OpenDrift framework is programmed in Python, and is based on four classes.
All classes are generic, and designed for specialisation by subclassing - ensuring flexibility, modularity and reusability.

* **Reader** takes care of retrieving driver/forcing data for requested positions/time from a given data source (model or observations)
* **LagrangianArray** describes a particle type (e.g. oil, larvae...) and its properties (position, mass, length...)
* **Writer** writes output to a specific file format at each step or at end of simulation

The main class (the "trajectory model") ties it all together and takes care of any other tasks:

* **OpenDriftSimulation**

  * Has an array of particles of a particular type (LagrangianArray)
  * Defines a spatial reference system (SRS), to which environmental data (from readers) will be interpolated
  * Methods for seeding/initialisation of particles.
  * Then for each time step (which might be negative):

    * Uses a set of Readers to retrieve driver/environmental data for the positions of the elements

      * includes reprojecting and interpolation of values from the Readers SRS to the common simulation SRS
    * Calls a purpose-specific method "update()" to update positions and properties of particles, based on environment data obtained by the readers.

      * Deactivates particles (stranded, evaporated etc).
    * Calls a Writer to append output to file.

  * Functionality to plot and analyse output/results (e.g as trajectories or concentration fields).

A subclass of OpenDriftSimulation defines (or imports/reuses) a particle type, and specifies a method to update particle properties based on environment variables (wind, waves, currents...). This function will focus only on the physical/chemical processes, and does not need to account for technical issues related to import, export, map projections etc. 

A class diagram and a more detailed description of the classes is given below.

.. image:: https://github.com/OpenDrift/opendrift/blob/master/docs/opendrift_UML.png?raw=true


Reader
*********

A "Reader" is a class to read/retrieve input/driver data from a particular source, or file format. This can be e.g. CF compliant NetCDF files, a specific model (ROMS) or source of observations (HF-radar, scatterometer...), with a specific file format/type, stored locally or retrieved through web protocols (e.g. OPeNDAP).
A Reader object must provide the following information, as reader by the class constructor (```__init__```):

* **proj4** (string)

  * declaring the spatial reference system ("projection") for which data can be requested
* **xmin, xmax, delta_x, ymin, ymax, delta_y** (floats)

  * defining grid bounds and resolution
  * x and y correspond to the given proj4 string
* **start_time, end_time, time_step** (datetime, timedelta)

  * defining temporal availability of variables
* **variables** (list of strings)

  * names of the variables which can be returned by this reader. Variable names must match names of the processes/functions of Model class (defined below) to allow generic use. This is ensured by using CF standard_name (whenever possible).

* **function get_variables(variables, time, x, y, depths, block)**

  * input:

    * variables is a list of variable names (CF standard_name) to be read
    * x, y, depth are arrays with coordinates of the positions for which environment data should be obtained.
    * block is a boolean, meaning:

      * True: a 3D block is returned for each requested variable, covering the requested x,y,depth space
      * False: variables returned are valid for the exact requested positions. If requested positions are not exactly matching the spatial grid (defined by min/max/delta x,y,depth) of the Reader, the nearest points are returned.

  * output:

    * a dictionary with:

      * key: variable name
      * value: 3D array (**depth, y, x**), or a scalar/constant. A list of values (3D or scalar) is returned for an ensemble model. Missing values (outside of domain, on land etc should be represented as masked values (NaNs or masked array).

    * arrays with corresponding values of depth, y, x


A few Readers are contained within the core OpenDrift repository, but (expert) users may make their own Readers for (local) data sources of interest. Therefore Readers class is defined to take care of only the very specific task of retrieving requested variables for requested position and time. More complex tasks such as reprojection and interpolation in time and space is taken care of centrally, preventing duplication of code.

LagrangianArray
******************

A "LagrangianArray" is the base class for any types of drifting elements. A LagrangianArray has a number of variables stored in an OrderedDict, with properties such as **dtype, unit, standard_name** and **default value**. The base class has only variables **lon, lat, depth**
Subclasses of LagrangianArray will add more properties, such as e.g. density and viscosity for oil elements, and mass and size (length/radius) for fish eggs, larvae or rigid objects. All properties of a LagrangianElements are arrays, representing a number of elements. Properties may also be 1D arrays, meaning that the same value represents all particles.

Model (subclass of OpenDriftSimulation)
*****************************************
A Model within the OpenDrift framework corresponds to the physics/chemistry/biology of a trajectory model in the classical sense (e.g. OD3D, Leeway, Ladim...). Any non-physical (technical) tasks are "outsourced" to other classes (Readers and Writers) and inherited/reused from classes LagrangianArray and OpenDriftSimulation.
A Model takes care of the processes affecting the particles (LagrangianArray) at each time step, such as:

* moving particles (x,y and vertical direction) in response to some forcing (e.g. current, wind, waves)
* changing any other properties of particles, as response to other physical proceses (evaporation, growth, capsizing...)
* check if particles are hitting land, and perform reponse (deactivating, rebouncing...)

All models must have a function update() which will be called each iteration of the trajectory simulation.
This method will use environment data (self.environment.x_sea_water_velocity, self.environment.y_sea_water_velocity, self.environment.x_wind, self.environment.sea_water_potential_temperature...) to update particle properties (self.elements.x, self.elements.depth, self.elements.mass...). The x-y coordinate frame is general, and the update function shall not "worry" about map projections and orientations etc.

Reader, LagrangianElement and OpenDriftSimulation are all `Abstract Base Classes <https://docs.python.org/2/library/abc.html>`_, meaning that they have to be subclassed/specialised to be used, and have to fulfil the mentioned capabilities.
