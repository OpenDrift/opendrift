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
     - no
     - no

   * - :mod:`OceanDrift3D <opendrift.models.oceandrift3D>`
     - yes
     - yes (optional wind_drift_factor)
     - yes
     - yes

   * - :mod:`OpenOil <opendrift.models.openoil>`
     - yes
     - yes (optional wind_drift_factor)
     - no
     - no

   * - :mod:`OpenOil3D <opendrift.models.openoil3D>`
     - yes
     - yes (optional wind_drift_factor)
     - yes
     - wave entrainment, turbulence, buoyancy

   * - :mod:`Leeway <opendrift.models.leeway>`
     - yes
     - yes, at an angle. category-based empirical empirical wind_drift_factor
     - implicit in wind
     - 

   * - :mod:`PelagicEgg <opendrift.models.pelagicegg>`
     - yes
     - no
     - yes
     - turbulence, buoyancy

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

Direct wind drift is only applied to elements/particles at the very surface. Elements may be seeded with a user defined property ```wind_drift_factor``` (default is typically 0.02, i.e. 2%) which determines the fraction of wind speed at which elements will be advected.

Surface Stokes Drift must be obtained from a wave model, whereas the depth dependency is parameterised according to `Breivik et al. (2014) <https://journals.ametsoc.org/doi/abs/10.1175/JPO-D-14-0020.1>`_.

Vertical entrainment, mixing and refloating is largely following `Tkalich and Chan (2003) <https://www.sciencedirect.com/science/article/pii/S0025326X02001789>`_ and `Röhrs et al. (2014) <https://onlinelibrary.wiley.com/doi/10.4319/lo.2014.59.4.1213/abstract>`_

All models are subclasses of `OpenDrift <opendrift.models.basemodel>`_, or `OpenDrift3D <opendrift.models.opendrift3D>`_ and inherits all core functionality from there. The OpenDrift class itself has no specification of advection or other processes, and can thus not be used directly.

As this page does not support wide tables, here is a separate table with additional comments:

.. list-table:: Title
   :widths: 20 40 40
   :header-rows: 1

   * - Model
     - Other processes
     - Purpose/comment

   * - OceanDrift
     - None
     - Basic model for e.g. drifters

   * - OpenOil
     - Oil weathering (evaporation, emulsification, biodegradation)
     - 2D oil model

   * - OpenOil3D
     - Oil weathering (evaporation, emulsification, biodegradation)
     - 3D oil model

   * - PelagicEgg
     - Biological behaviour (TBD)
     - Fish eggs

   * - OpenBerg
     - 
     - Ice bergs
