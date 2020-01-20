How to write a new module
==========================

The best way to get started in developing a new model, is to copy and modify an existing working model/module (e.g. `oceandrift <https://github.com/opendrift/opendrift/blob/master/opendrift/models/oceandrift.py>`_) or look at the `model_template <https://github.com/opendrift/opendrift/blob/master/opendrift/models/model_template.py>`_ which contains instructions/comments. A general explanation is given below.


The main class of the OpenDrift framework is the **OpenDriftSimulation** class in the module `basemodel.py <https://github.com/opendrift/opendrift/blob/master/opendrift/models/basemodel.py>`_, stored in the folder `models <https://github.com/opendrift/opendrift/blob/master/opendrift/models/>`_.

A purpose specific trajectory model is made by subclassing **OpenDriftSimulation** (or any of its subclasses/models), mainly by:

  * defining (or importing/reusing) a particle/element type (e.g. oil particle, larvae particle...), by specifying names and default values of the (additional) properties the elements shall have.
  * defining which environment variables (wind, waves, current...) shall be needed by the model to update positions and properties of elements (below).

    * these environment variables are obtained by :doc:`Readers <data_model>`, however, when writing a model/module one do not need to worry about how environment variables are obtained.
  * defining a class method **update()** to be called at each model time step.

The method **update()** will in principle do three different things, regardless of the specific application:

1. Moving elements
##################

The positions of elements are stored as the element properties ``longitude`` and ``latitude`` with units of degrees, as well as the vertical coordinate ``z`` with units of meters, negative in the ocean. These properties might in principle be modified directly within the update() method, however, it is both safer and simpler to rather use the built in convenience methods, such as::

    self.update_positions(x_velocity, y_velocity)

This method takes care of e.g. rotation of vectors from the specific x-y coordinate system of the simulation (Pyproj projection) to determine the correct update of the position (longitude and latitude), also taking the simulation time_step into account.
The ``x_velocity`` and ``y_velocity`` might be any velocity provided by the readers (currents, wind, stokes drift...) possibly scaled by the model application (e.g. 0.02 times the wind, which is often used for wind drift at ocean surface).

By inheritance, models may also reuse higher level advection methods inherited from their superclasses. E.g. the model `OceanDrift <https://github.com/opendrift/opendrift/blob/master/opendrift/models/oceandrift.py>`_) uses very simple methods which takes care of everything needed for a correct and efficient advection::

    # Simply move particles with ambient current
    self.advect_ocean_current()

    # Advect particles due to wind drag (according to specified wind_drift_factor)
    self.advect_wind()

Notice that these (convenience) methods are inherited in the specific model `OceanDrift <https://github.com/opendrift/opendrift/blob/master/opendrift/models/oceandrift.py>`_ from the generic model `OpenDrift <https://github.com/opendrift/opendrift/blob/master/opendrift/models/opendrift.py>`_, which again inherits them from the helper class `PhysicsMethods <https://github.com/opendrift/opendrift/blob/master/opendrift/models/physics_methods.py>`_).
Reusing functionality by inheritance makes it simpler and faster to develop and maintain new models and functionality.
The vertical coordinate (z) is modified e.g. by common turbulence schemes.

2. Modifying element properties
###############################

Element properties may be modified freely and directly within the method **update()**. The main "resources" available are the Python recarrays:

 * ``self.elements`` which has any element properties (mass, length, orientation, age...) as attributes, and
 * ``self.environment`` which has any environment variables (x_wind, y_wind, sea_water_temperature) as attributes.

E.g. if one (for some reason) would like to set the length of the elements equal to half of the seawater temperature::

    self.elements.length = self.environment.sea_water_temperature / 2

Inheritance is also used for modifying properties, e.g. `OpenOil3D <https://github.com/opendrift/opendrift/blob/master/opendrift/models/openoil3D.py>`_ which inherits the full oil weathering from `OpenOil <https://github.com/opendrift/opendrift/blob/master/opendrift/models/openoil.py>`_, and simply calls the method ``self.oil_weathering()`` within its own ``update()`` method (in addition to the vertical advection not found in the 3D version OpenOil).

3. Deactivating elements
########################

Deactivation of elements must always be done with the special class method ``deactivate_elements(indices, reason)``
Elements may be deactivated for any reason as determined by the model developer, and an arbitrary text string is provided as explanation. E.g. if one want to deactivate all elements whose mass is below, say 0.5 kg::

    self.deactivate_elements(self.elements.mass < 0.5, reason='vanished')

``self.elements.mass`` is a numpy array with one value per element/particle, and thus ``self.elements.mass < 0.5`` returns the indices of the elements for which the mass is below 0.5 kg (units are specified in element type declaration above).
The text string ``reason`` serves as an explanation of why the element was deactivated, and is used e.g. for labeling figures and when analysing results. The element property ``status`` will be updated accordingly, as the element is deactivated. All elements have status 0 (integer) as default, meaning ``active``. Each time (in method **update()**) a new deactivation reason is used, the string (reason) is added to the list ``status_categories``. The list of deactivation reason is also added to the metadata in exported netCDF files, e.g::

    int status(trajectory, time) ;
		    status:flag_values = 0, 1, 2 ;
		    status:flag_meanings = "active evaporated stranded" ;

where one can see that elements have been deactivated for two different reasons ('evaporated' **1**, and 'stranded' **2**).
