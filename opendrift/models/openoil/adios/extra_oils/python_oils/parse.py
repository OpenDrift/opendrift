import NO00159_SVALE as d
import adios_db.scripting as ads
from adios_db.models.oil.values import Reference
from adios_db.models.oil.physical_properties import DynamicViscosityPoint
from adios_db.models.oil.sara import Sara
from adios_db.models.oil.distillation import DistCutList


oil = ads.Oil(d.id)
oil.metadata.name = d.name
oil.metadata.API = (141.5/d.specific_gravity) - 131.5
oil.metadata.product_type = d.product_type
oil.metadata.reference = Reference(year=d.reference_year, reference=d.reference)

oil.sub_samples.extend([ads.Sample(), ads.Sample(), ads.Sample(), ads.Sample()])
for i, (ss, w, vol, visc, dens) in enumerate(
        zip(oil.sub_samples, d.weathering_names, d.weathering_volume,
            d.weathering_visc_mpas, d.weathering_density)):
    print(i, w, vol, visc, dens)
    ss.metadata.name = w
    ss.metadata.short_name = w
    ss.metadata.fraction_evaporated = ads.MassFraction(100-vol, unit='%')
    ss.physical_properties.densities.append(
            ads.DensityPoint(density = ads.Density(value=dens, unit="kg/m^3"),
            ref_temp=ads.Temperature(value=d.weathering_ref_temp, unit='C')))
    ###############################
    # TODO check mpas -> kg/ms
    ###############################
    ss.physical_properties.dynamic_viscosities.append(DynamicViscosityPoint(
            viscosity=ads.DynamicViscosity(value=visc/100, unit='kg/(m s)'),
            ref_temp=ads.Temperature(value=d.weathering_ref_temp, unit='C')))

    print(ss)

for i, s in enumerate(d.weathering_names):
    oil.sub_samples[i].metadata.name = s
    oil.sub_samples[i].metadata.short_name = s

oil.sub_samples[0].SARA = Sara(asphaltenes=ads.MassFraction(value=d.asphaltenes_percent, unit='%'))
oil.sub_samples[0].bulk_composition.append(ads.MassFraction(value=d.wax_percent, unit='%'))

cut_temps = [c[0] for c in d.cuts]
cut_frac = [c[1] for c in d.cuts]
print(cut_temps)
print(cut_frac)
cutl = DistCutList.from_data_arrays(temps=cut_temps, temp_unit='C',
                                    fractions=cut_frac, frac_unit='%')
oil.sub_samples[0].distillation_data.cuts = cutl
oil.sub_samples[0].distillation_data.type = 'volume fraction'
oil.sub_samples[0].distillation_data.fraction_recovered = ads.Concentration(max_value=1.0, unit="fraction")


print(oil)
oil.status = oil.validate()
print(oil.py_json())
print('\n')
print(oil.sub_samples[0])
print('\n')

for s in oil.status:
    print(s)

