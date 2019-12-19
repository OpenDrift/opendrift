Graphical User Interface
========================

OpenDrift comes with a basic Graphical User Interface, which can be started with the command:

.. code::

   $ opendrift_gui.py

given that the `scripts-folder <https://github.com/OpenDrift/opendrift/tree/master/opendrift/scripts>`_ is added to your PATH-variable. Otherwise, the command `$ python opendrift_gui.py`
should be run from within the scripts-folder. On Windows, the opendrift_gui.py script may be launched by a double-click.

.. image:: images/opendrift_gui.png

Choosing a model and object/oiltype
-----------------------------------
Presently, only the Leeway model (Search and Rescue) and OpenOil (oil drift) are available in the GUI.

The first step is to select either Leeway or OpenOil from the upper menu.
Based on this choice, the second dropdown-meny will allow to choose either
  * an oiltype (from `NOAA OilLibrary <https://github.com/NOAA-ORR-ERD/OilLibrary>`_), or
  * a `Leeway object category <https://www.dropbox.com/s/d5jhwo6y5zldxgh/Drivbanehjelp.pdf?raw=1>`_.

The following steps are the same for both models:

Seeding of particles
--------------------
* Choose a start point for the seeding (initialisation) of particles. Latitude and longitude must be given in decimal degrees, and negative longitudes correspond to west of Greenwich.
* The radius is an uncertainty radius (standard deviation, given in meters) around the given location within which all the (5000) elements will be seeded.
* Optional:
  * the end time can be changed if a continuous release of particles over a time span is desired, simulating e.g. a leaking well,
  * the end position can be changed if particles shall be released between the start and end positions, simulating e.g. a ship track.
  * or a combination of both.
* Note that whenever the start position or time is updated, the end position/time will be copied/overwritten.
* By clicking `Check seeding` it is possible to see a plot of the seedig positions before the simulation is started. Note that it takes about 5 seconds to generate this plot, and even 30 seconds if full resolution coastline is chosen.

Simulation time
---------------
* Simulations may by run forwards or backwards (backtracking, to check from where the object or oil has drifted)
* The simulation duration is given in hours. Note that simulations longer than 100 hours may take significant time (> 10 minutes), especially for oil drift.

Forcing data
------------
The following sources are used for current, wind and waves, in prioritised order:

* **Current**:

  - MET Norway `NorKyst ROMS model <http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be>`_ covering the full Norwegian coastline at 800 m resolution. Coverage from 14 Sept 2016 until about 60 hours into the future
  - MET Norway `Nordic4 ROMS model <http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be>`_ covering the North Atlantic 4 km resolution. Coverage from 19 Sept 2016 until about 5 days into the future
  - `Global FOAM model of UK Met Office <http://data.ncof.co.uk/thredds/dodsC/METOFFICE-GLO-AF-PHYS-HOURLY-CUR>`_ at 25 km resolution. Coverage from 11 July 2017 until 12 days into the future.

* **Wind**:

  - MET Norway `MEPS <http://thredds.met.no/thredds/dodsC/meps25files/meps_det_pp_2_5km_latest.nc>`_ (Arome) model covering Norway and Sweden with surrounding ocean at 2.5 km resolution. Coverage from a few hours back until about 60 hours into the future.
  - MET Norway `Arome Arctic <http://thredds.met.no/thredds/dodsC/aromearcticlatest/arome_arctic_pp_2_5km_latest.nc>`_ model covering the Barents Sea and part of the Arctic at 2.5 km resolution. Coverage from a few hours back until about 60 hours into the future.
  - `Global NCEP GFS <http://oos.soest.hawaii.edu/thredds/dodsC/hioos/model/atm/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd>`_ winds. Coverage from 6 May 2011 until 7 days into the future.

* **Waves**:

  - MET Norway `MyWaveWAM model <http://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be>`_ at 8 km resolution, covering the North Atlantic. Coverage from 18 Sept 2016 until about 60 hours into the future.

Whenever wave data is not available, OpenDrift will calculate wave parameters from wind.
The forcing data are obtained/downloaded from remote servers (Thredds, OPeNDAP) during the simulation.
Possibility to modify/configure the forcing data is planned for the near future.

Simulation and analysis
-------------------------
* The simulation is started by clicking `PEIS PAO`, which is vossamaol for `START`.

* The simulation will take from less than a minute for a short (12 hours) Leeway simulation and up to 10 minutes or more for a long (100 hours) oil drift simulation. Oil drift simulations takes more time due to heavy vertical turbulence calculations.

  * A progress-bar is not yet available in the GUI, but is planned for the near future.
  * It is presently not possible to stop an ongoing simulation, except for using a task manager (or `kill` command in Linux). The user interface is also unresponsive during the simulation.

* After the simulation, the result may be analysed with buttons appearing at the bottom:

  * `Animation`: an animation of the simulation is shown. There is presently no possibility to save this animation to file from the GUI.
  * `Oil budget` (oil drift only): the time-history of the :mod:`oil budget <opendrift.models.openoil>` is shown
  * MET Norway only: a button `Show in Diana` will launch the Diana analysis tool.

    * Note that this is local Diana, not VGL-version.
    * The present and previous simulations may be viewed by opening `Felter -> Drivbane_OpenDrift` and then selecting `Leeway` or `Oljedrift`. Note that the simulations will be sorted by the starting-time of the simulations, not the time at which the simulation was performed. Thus if you performed a simulation starting back in time, it might not be the last item in the list (`Referansetid`). Click on all the items in felter (submerged, stranded, surface), and then `Utfør`. For oil drift simulations, submerged particles are blue, surface oil is grey or black.

