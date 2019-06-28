API Reference for Opendrift config file
==============================================

The opendrift configuration file consists of a yaml file including details on 
 - the simulation to run
 - the "scheduler" setup which defines the meta parameters related to when to run, time limits, which docker to use or how many cores to allocate etc..

This documents focuses on the parameters used to setup the OpenDrift simulation. The scheduler part is described in full details here: https://github.com/metocean/scheduler/blob/master/docs/actions.rst

a wide range of parameters can be  passed to the opendrift object's methods using kwargs via :
extra_model_args
extra_seed_args
readers' options
extra_output_args

**Opendrift Model**

::

    imp: test_OpenOil3D_CMEMS # implementation name
    model: OpenOil3D # opendrift object to use in simulation
    extra_model_args: {logfile: 'test_OpenOil3D_CMEMS.log' , weathering_model: 'noaa'}
    rootdir: 'C:\github\opendrift\examples_msl' # where the simulation will be run

**Release timing, position, and geographic frame**

::

    # release
    nb_parts: 1000
    start_time: '01-01-2017 00:00'
    end_time:   '08-01-2017 00:00'
    # duration_hours: 24.0 #= int(self.durationhours.get())*3600/time_step    
    end_time_run: '08-01-2017 00:00' #not used if duration is specified
    time_step_sec: 900 #= 900  # 15 minutes
    time_step_sec_output: 1800 #= timedelta(minutes=30)
    position: {lon: 174.1656,lat: -40.3346, radius: 10.,z: 0.0} 
    # z can be 'seafloor' or 'seafloor+2' or 'surface' , or scalar
    # optional; final lon,lat,z specified as
    # end_position: {elon: 174.1656,elat: -40.3346, eradius: 10.,ez: 0.0} 
    extra_seed_args: {} # these will be passed as **kwargs to seed_elements() opendrift's method
    # e.g. extra_seed_args: {objectType: 26,} # 26 = Life-raft, no ballast
    # extra_seed_args: {terminal_velocity: -0.001, wind_drift_factor: 0.02}
    model_frame:    {llcrnrlon: 168.0199, llcrnrlat: -42.8449 , urcrnrlon: 177.4601 , urcrnrlat: -37.9051}
    # geographical extents of simulation
    basemap_resolution: 'h' 
    # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)

**Readers (i.e. hydrodynamic/wind/wave fields)**

This block specifies the metocean (hydrodynamic/wind/wave) data to use for simulations. Available options for now include specifying path to a local file, download data from Copernius Marine Service, internal MetOcean UDS, point to a THREDDS server (Opendrift functionality).
Note the order in which readers are specified will define the reader's priority. In example below, model will try to get data first from ocean0, then ocean1, then ocean2. Same can be reproduced for wind and wave readers.

::
    
    readers:
        ocean0: # reader name used internally in opendrift
            reader_type: [reader_netCDF_CF_generic] #reader_type can be any found in 'from opendrift.readers'  
            cmems_download: True # data needs to be downloaded from Copernicus Marine Service
            dset:    ['GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'] # dataset to use
            vars:    ['utotal','vtotal'] # variables that will be queried
            datatype: 'hc' # optional 
            timestep: 1.0 # not used for CMEMS queries
            boundary:    [170.0,176.0,-42.0,-38.0] 
            # boundary is optional - if not input use a square of 2x2 deg around release
            options: {use_log_profile: True, z0: 0.001}  # options to be passed as kwargs to reader's __init__() function - optional
        ocean1: # 
            reader_type: [reader_netCDF_MetOcean] #reader_type can be any found in 'from opendrift.readers'  
            # filename: ['F:/metocean/R&D_OpenDrift/benchmark_runs/cnz20040101_00z_surf.nc'] # can be many files
            udshost: 'http://uds1.rag.metocean.co.nz:9191/uds' # data needs to be downloaded from UDS
            dset:    ['roms_cnz_surf'] # dataset to use
            vars:    ['um','vm','dep'] # variables to query
            datatype: 'hc' # datatype
            timestep: 1.0 # required timestep
            boundary:    [173.0,175.5,-42.0,-39.0] # optional - if not input use a squar of 2x2 deg.
            options: {use_log_profile: True, z0: 0.001} # options to be passed as kwargs to reader's __init__() function - optional
        ocean2: # 
            reader_type: [reader_netCDF_MetOcean] #reader_type can be any found in 'from opendrift.readers'  
            filename: ['F:/metocean/R&D_OpenDrift/benchmark_runs/cnz20040101_00z_surf.nc'] # path to local file (can use wildcards)
        wave: 
            reader_type: [reader_netCDF_CF_generic] #reader_type can be any found in 'from opendrift.readers'  
            cmems_download: True # data needs to be downloaded from UDS
            dset:    ['GLOBAL_ANALYSIS_FORECAST_WAV_001_027-TDS']
            vars:    ['VHM0','VTPK','VTM10','VPED','VMDR','VSDX ','VSDY']
            datatype: 'hc'
            timestep: 1.0
            boundary:    [170.0,176.0,-42.0,-38.0] # optional - if not input use a squar of 2x2 deg.
        wind: 
            reader_type: [reader_netCDF_CF_generic] #reader_type can be any found in 'from opendrift.readers'  
            cmems_download: True # data needs to be downloaded from UDS
            dset:    ['WIND_GLO_WIND_L4_REP_OBSERVATIONS_012_006-TDS'] #WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004
            vars:    ['eastward_wind','northward_wind']
            datatype: 'hc'
            timestep: 1.0
            boundary:    [170.0,176.0,-42.0,-38.0] # optional - if not input use a squar of 2x2 deg.
        bathy: # workaround to access depth information even when using CMEMS data - query timeseries of gridded depth from UDS
            reader_type: [reader_netCDF_MetOcean] #reader_type can be any found in 'from opendrift.readers'  
            udshost: http://uds1.rag.metocean.co.nz:9191/uds # data needs to be downloaded from UDS
            dset:    ['roms_nz_surf']
            vars:    ['dep']
            datatype: 'hc'
            timestep: 1.0
            boundary:    [170.0,176.0,-42.0,-38.0] # optional - if not input use a squar of 2x2 deg.
            # options: { always_valid: True } # will be passed as kwargs to reader's __init__()

   fallback_values: {ocean_vertical_diffusivity: 0.0001, x_wind: 0.0, y_wind: 0.0} 
   # used to define constant/fall back values in case no data can be interpolated from readers


**Opendrift model configuration**


All details of the Opendrift model configuration can be set here. They will overwrite default values. All available configuration parameters are given in opendrift source code, at begining of each ./models/model.py scripts

::

    # if nothing specfied, default config is used 
    config:
        general:
            coastline_action: 'stranding' # option('none', 'stranding', 'previous', default='stranding')
        seed:
            ocean_only: True #boolean(default=True)
            # oiltype: 'GULLFAKS, EXXON' - can be set here instead of extra_seed_args 'ARABIAN MEDIUM, API' 
        drift:
            scheme: 'runge-kutta4' #option('euler', 'runge-kutta4', 'runge-kutta4', default='euler')
            current_uncertainty: 0.025
            wind_uncertainty: 2.0
            # max_age_seconds: float(min=0, default=None)
            stokes_drift: True #boolean(default=True)
            # wind_drift_depth: 3.0
            # current_uncertainty: float(min=0, max=5, default=0)
            # current_uncertainty_uniform: float(min=0, max=5, default=0)
            # wind_uncertainty: float(min=0, max=5, default=0)
            # relative_wind: boolean(default=False)
            # lift_to_seafloor: boolean(default=True)
            # truncate_ocean_model_below_m: float(min=0, max=10000, default=None)
            # deactivate_north_of: float(min=-90, max=90, default=None)
            # deactivate_south_of: float(min=-90, max=90, default=None)
            # deactivate_east_of: float(min=-360, max=360, default=None)
            # deactivate_west_of: float(min=-360, max=360, default=None)
            # use_tabularised_stokes_drift: boolean(default=False)
            # tabularised_stokes_drift_fetch: option(5000, 25000, 50000, default=25000)
        processes:
            turbulentmixing: True
            verticaladvection: False
            dispersion: True
            evaporation: True
            emulsification: True
            update_oilfilm_thickness: False
        # input:
        #     [[spill]]
        #         oil_type = option(%s, default=%s)
        #         droplet_diameter_min_subsea = float(min=1e-8, max=1, default=0.0005)
        #         droplet_diameter_max_subsea = float(min=1e-8, max=1, default=0.005)
        wave_entrainment:
            droplet_size_distribution: 'Johansen et al. (2015)' #= option('Exponential', 'Johansen et al. (2015)', 'Li et al. (2017)', default='Johansen et al. (2015)')
            entrainment_rate: 'Li et al. (2017)' #= option('Tkalich & Chan (2002)', 'Li et al. (2017)', default='Li et al. (2017)')
        turbulentmixing:
            timestep: 5 #= float(min=0.1, max=3600, default=60.)
            verticalresolution: 1. #= float(min=0.01, max=10, default = 1.)
            diffusivitymodel: 'environment' #= option('environment', 'stepfunction', 'windspeed_Sundby1983', 'windspeed_Large1994', 'gls_tke', default='environment')
            TSprofiles: False #= boolean(default=False)
            droplet_diameter_min_wavebreaking: 1e-5 #= float(default=1e-5, min=1e-8, max=1)
            droplet_diameter_max_wavebreaking: 2e-3 #= float(default=2e-3, min=1e-8, max=1)
            droplet_size_exponent: 0 #= float(default=0, min=-10, max=10)

**Run details**

::

    run_backwards: False
    stop_on_error: False
    outfile: output_cmems_openoil3d.nc #self.outputdir + '/opendrift_'  # could be defined by default from implemtation name ?
    extra_output_args: {} # additional outputs arguments to be passed as kwargs to opendrift run() method

**Post Processing**

::

    post_process:
        show_anim:  False
        save_anim:  False
        show_plot:  True
        save_plot:  False
        show_oil_budget: True
        save_oil_budget: False
        # more to come


**Notes**
--------------------------
Any parameters can be  passed to the opendrift object's methods using kwargs via

::    

   extra_model_args # passed to __init_()

::    

   extra_seed_args # passed to seed_elements()

::    

   options (for readers) # passed to Reader.__init_()

::    

   extra_output_args # passed to run()


see examples full configuration files in https://github.com/simonweppe/opendrift/tree/master/examples_msl/config
