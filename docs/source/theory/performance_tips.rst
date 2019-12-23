Improving performance
=====================

OpenDrift is fairly optimised Python code, but simulations may take long time (even several days) if many elements (~millions) and many time steps (~thousands).
We hope to parallelise OpenDrift in the near future, to improve performance on machines with many cores. In the meantime, the following can be done to increase performance or decrease memory and disk usage:

* Decrease the number of particles to a minimum
* Not using as larger time steps than necessary. When calling method ``run()``, one can specify ``time_step`` which is the internal time step for calculation, and ``time_step_output`` which is the steps at which data are stored in memory (and written to file if ``outfile`` is specified)
* Reduce memory consumption by either:

  * specifying ``export_buffer_length`` in run(). 100 is default, but a smaller value may be specified if there are very many elements (~millions)
  * limit the number of variables saved to disk (and stored in memory) by providing a list of needed variables as the list ``export_variables`` when calling run()

* Reading data from local files is faster than from remote servers, e.g. Thredds, although the latter is very convenient.
* For the 3D models, vertical mixing takes considerable time. If not needed, this can be switched off by setting  ``o.config['processes']['turbulentmixing'] = False``

  * An upper limit may be given to the number of iterations per outer timestep, as illustrated in example_long_verticalmixing.py.
  * The mixing can also be significantly speeded up by using constant TS-profiles during the iterations, by specifying ``o.config['turbulentmixing']['TSprofiles'] = False`` (which is default)
  * The mixing may also be speeded up by setting a larger vertical resolution, which further allows for a larger time step (order dz^2).

    * mixing vertical resolution can be adjusted by ``o.config['turbulentmixing']['verticalresolution'] = 2 # meters``
    * mixing time step can be adjusted by ``o.config['turbulentmixing']['timestep'] = 4  # seconds``
    * The time step should be small enough to avoid debug-warnings like ``WARNING! some elements have p+q>1`` A vertical resolution of 2 m seems to allow for a time step of ~4 seconds.
* Quasi-parallelisation can be done by splitting the initial elements up in geographical sub-regions, and running simulations in parallel for each subset.


Some notes about performance
*******************************

Each OpenDrift module has an attribute ``max_speed`` which indicates the maximum velocity that elements will likely encounter. This is accessible as ``o.max_speed`` for an OpenDrift object ``o``. Typical values are 2 m/s for ocean drift models, but as large as 12 m/s for the WindBlow module where elements move at the speed of the wind. It might be overriden by the user at any time, simply by writing ``o.max_speed = 3`` before starting a run. This parameter is not an absolute bound for instantaneous velocity, but rather for the average velocity in a particular direction over a time period. It serves two purposes:

* If a reader providing landmask (``land_binary_mask``) has not been added to an OpenDrift object before a run is started, a GSHHG vector landmask is created based on an estimated bound of where elements might move. This landmask will cover a rectangle around the positions at which elements have been seeded, with a spatial buffer to east/west/north/south equal to ``o.max_speed
  * simulation_duration_seconds``. It is of interest to keep this buffer as small as possible, as checking land encounter for many particles versus a large area of coastline and islands may be costly.

* During the run, all readers are called to return data (2D or 3D "blocks") covering all the elements. This block of data should not only cover elements at their present position, but also the area over which elements might move within the next timestep of the reader (e.g. typically 1, 3, 6 or 24 hours for an ocean model). As with the landmask above, this buffer around present positions is calculated as ``o.max_speed * timestep_of_reader``. Again, it of interest to keep the buffer as small as possible to save time reading data from file (or from Thredds), but still avoiding that elements will leave the reader-block coverage within the timestep.  
