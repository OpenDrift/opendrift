How to choose which model to use
================================

OpenDrift contains a few specific models for several applications such as oil drift, search and rescue and fish eggs.

The table below shows an overview of the advection processes within the main models:

Model | Move with currents | Direct wind drift | Stokes drift  |  Vertical motion |
--- | --- | --- | --- | --- |
[**OceanDrift**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/oceandrift.py) | yes | yes (optional wind_drift_factor) | no | no |
[**OceanDrift3D**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/oceandrift3D.py) | yes | yes (optional wind_drift_factor) | yes | yes |
[**OpenOil**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/openoil.py) | yes | yes (optional wind_drift_factor) | no | no |
[**OpenOil3D**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/openoil3D.py) | yes | yes (optional wind_drift_factor) | yes | wave entrainment, turbulence, buoyancy|
[**Leeway**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/leeway.py) | yes | yes, at an angle. category-based empirical empirical wind_drift_factor | no (included in wind drift) | no |
[**PelagicEgg**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/pelagicegg.py) | yes | no | yes | turbulence, buoyancy |
[**PlastDrift**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/plastdrift.py) | yes | no | yes | empirical-statistical depth, exponential decrease with depth, depending on turbulence/wind |

Direct wind drift is only applied to elements/particles at the very surface. Elements may be seeded with a user defined property ```wind_drift_factor``` (default is typically 0.02, i.e. 2%) which determines the fraction of wind speed at which elements will be advected.

Surface Stokes Drift must be obtained from a wave model, whereas the depth dependency is parameterised according to [Breivik et al. (2014)](http://journals.ametsoc.org/doi/abs/10.1175/JPO-D-14-0020.1).

Vertical entrainment, mixing and refloating is largely following [Tkalich and Chan (2003)](http://www.sciencedirect.com/science/article/pii/S0025326X02001789) and [Röhrs et al. (2014)](http://onlinelibrary.wiley.com/doi/10.4319/lo.2014.59.4.1213/abstract)

All models are subclasses of [**OpenDrift**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/basemodel.py), or [**OpenDrift3D**](https://github.com/knutfrode/opendrift/blob/master/opendrift/models/opendrift3D.py) and inherits all core functionality from there. The OpenDrift class itself has no specification of advection or other processes, and can thus not be used directly.

***
As GitHub wiki does not support wide tables, here is a separate table with additional comments:

Model | Other processes  | Purpose/comment
--- | --- | --- |
**OceanDrift** | None | Basic model for e.g. drifters|
**OpenOil** | Oil weathering (evaporation and emulsification, but presently very simple, to be improved) | 2D oil model|
**OpenOil3D** | Oil weathering (evaporation and emulsification, but presently very simple, to be improved) | 3D oil model|
**Leeway** | Jibing (shifting orientation) 4% probability per hour | Search and rescue |
**PelagicEgg** | Biological behaviour (TBD)| Fish eggs |
