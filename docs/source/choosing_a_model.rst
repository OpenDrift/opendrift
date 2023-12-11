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

   * - :mod:`OpenBergOld <opendrift.models.openberg_old>`
     - yes
     - yes
     - implicit in wind
     - no

Direct wind drift is only applied to elements/particles at the very surface. Elements may be seeded with a user defined property ``wind_drift_factor`` (default is typically 0.02, i.e. 2%) which determines the fraction of wind speed at which elements will be advected.

Surface Stokes Drift must be obtained from a wave model, whereas the depth dependency is parameterised according to `Breivik et al. (2014) <https://journals.ametsoc.org/doi/abs/10.1175/JPO-D-14-0020.1>`_.

Vertical entrainment, mixing and refloating is largely following `Röhrs et al. (2018) <https://doi.org/10.5194/os-14-1581-2018>`_

All models are subclasses of :py:mod:`OpenDrift <opendrift.models.basemodel>`, or :py:mod:`OceanDrift <opendrift.models.oceandrift>` and inherits all core functionality from there. The OpenDrift class itself has no specification of advection or other processes, and can thus not be used directly.
