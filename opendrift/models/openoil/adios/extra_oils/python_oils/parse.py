def parse():
    import NO00160_BREIDABLIKK as d
    import NO00159_SVALE as d
    import NO00161_STAER as d
    import NO00162_YME as d
    import NO00163_OFELIA as d
    import NO00164_HEIDRUN_AARE as d
    import adios_db.scripting as ads
    from adios_db.models.oil.values import Reference
    from adios_db.models.common.measurement import (MassFraction,
                                                    Temperature,
                                                    MassOrVolumeFraction,
                                                    AnyUnit)

    from adios_db.models.oil.physical_properties import DynamicViscosityPoint, PourPoint, FlashPoint
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
        # TODO check mpas -> kg/ms -> checked, 1000 should be correct
        ###############################
        ss.physical_properties.dynamic_viscosities.append(DynamicViscosityPoint(
                viscosity=ads.DynamicViscosity(value=visc/1000, unit='kg/(m s)'),
                ref_temp=ads.Temperature(value=d.weathering_ref_temp, unit='C')))
        if d.pour_point is not None:
            ss.physical_properties.pour_point = PourPoint(measurement=Temperature(value=d.pour_point, unit='C'))
        if d.flash_point is not None:
            ss.physical_properties.flash_point = FlashPoint(measurement=Temperature(value=d.flash_point, unit='C'))

    for i, s in enumerate(d.weathering_names):
        oil.sub_samples[i].metadata.name = s
        oil.sub_samples[i].metadata.short_name = s

    oil.sub_samples[0].SARA = Sara(asphaltenes=ads.MassFraction(value=d.asphaltenes_percent, unit='%'))
    oil.sub_samples[0].bulk_composition.append(ads.MassFraction(value=d.wax_percent, unit='%'))

    cut_temps = [c[0] for c in d.cuts]
    cut_frac = [c[1] for c in d.cuts]
    cutl = DistCutList.from_data_arrays(temps=cut_temps, temp_unit='C',
                                        fractions=cut_frac, frac_unit='%')
    oil.sub_samples[0].distillation_data.cuts = cutl
    oil.sub_samples[0].distillation_data.type = 'volume fraction'
    oil.sub_samples[0].distillation_data.fraction_recovered = ads.Concentration(value=1.0, unit="fraction")

    oil.status = oil.validate()

    for s in oil.status:
        print(s)

    print(oil.status)
    print(oil.metadata.gnome_suitable, 'GNOME suitable')
    oil.to_file(f'{d.id}.json')


if __name__ == '__main__':
    parse()

