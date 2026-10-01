How to choose which model to use
================================

OpenDrift contains a few specific models for several applications such as oil drift, search and rescue and fish eggs.

The table below shows an overview of the advection processes within the main models:

.. list-table::
   :widths: 20 10 30 20 20
   :header-rows: 1

   * - Model
     - Move with currents
     - Direct wind drift
     - Stokes drift
     - Vertical motion

   * - :mod:`OceanDrift <opendrift.models.oceandrift>`
     - yes
     - yes (optional wind_drift_factor)
     - yes
     - advection, turbulence

   * - :mod:`OpenOil <opendrift.models.openoil>`
     - yes
     - yes (optional wind_drift_factor)
     - yes
     - wave entrainment, turbulence, buoyancy, advection

   * - :mod:`Leeway <opendrift.models.leeway>`
     - yes
     - yes, at an angle. category-based empirical empirical wind_drift_factor
     - implicit in wind
     - no

   * - :mod:`PelagicEgg <opendrift.models.pelagicegg>`
     - yes
     - no
     - yes
     - turbulence, buoyancy, advection

   * - :mod:`PlastDrift <opendrift.models.plastdrift>`
     - yes
     - no
     - yes
     - empirical-statistical depth, exponential decrease with depth, depending on turbulence/wind

   * - :mod:`OpenBerg <opendrift.models.openberg>`
     - yes
     - yes
     - implicit in wind
     - no

Direct wind drift is only applied to elements/particles at the very surface. Elements may be seeded with a user defined property ``wind_drift_factor`` (default is typically 0.02, i.e. 2%) which determines the fraction of wind speed at which elements will be advected.

By default the wind drift is in the same direction as the wind. The configuration option ``drift:wind_drift_angle`` (degrees, default 0) rotates the wind before the wind drift is calculated. A positive angle turns the drift to the right of the wind in the Northern Hemisphere and to the left in the Southern Hemisphere, due to the Coriolis effect. The rotation doesn't change the wind speed. Only a small angle should be used, since the wind-driven turning of the surface current (Ekman drift) is normally already included in the ocean model current (`Jones et al. (2016) <https://doi.org/10.1002/2016JC012113>`_), and a large angle would apparently count this turning twice. A small angle may still be useful, because the uppermost layer of an ocean model can underestimate the current at the very surface (Jones et al., 2016). For example, `Reed et al. (1994) <https://doi.org/10.1016/1353-2561(94)90009-4>`_ found a best fit to an experimental oil spill with a deflection angle of about 3 degrees to the right.

Surface Stokes Drift must be obtained from a wave model, whereas the depth dependency is parameterised according to `Breivik et al. (2014) <https://journals.ametsoc.org/doi/abs/10.1175/JPO-D-14-0020.1>`_.

Vertical entrainment, mixing and refloating is largely following `Röhrs et al. (2018) <https://doi.org/10.5194/os-14-1581-2018>`_

All models are subclasses of :py:mod:`OpenDrift <opendrift.models.basemodel>`, or :py:mod:`OceanDrift <opendrift.models.oceandrift>` and inherits all core functionality from there. The OpenDrift class itself has no specification of advection or other processes, and can thus not be used directly.
