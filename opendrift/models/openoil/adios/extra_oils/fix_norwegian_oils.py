import glob
from datetime import datetime
import numpy as np
import argparse
import adios_db
from adios_db.scripting import Concentration
from adios_db.models.oil.metadata import MetaData, ChangeLogEntry, ChangeLog
from adios_db.models.oil.oil import Oil
from adios_db.models.oil.physical_properties import DynamicViscosityList

def fix_norwegian_oils():
    parser = argparse.ArgumentParser(description='Fix JSON oil files')
    parser.add_argument('-i', metavar='--input', type=str, help='Input folder with json files', default='none')
    parser.add_argument('-o', metavar='--output', type=str, help='Output folder to save updated json files', default='none')
    args = parser.parse_args()
    if args.i == 'none':
        parser.print_help()

    folder_in = args.i + '/'
    folder_out = args.o + '/'

    print(folder_in, folder_out)

    prefix = 'NO*.json'
    oils = glob.glob(folder_in + prefix)

    for oil in oils:
        o = Oil.from_file(oil)
        eb = len(o.validate())

        cl = ChangeLogEntry(name='Knut-Frode Dagestad', date=datetime.now().isoformat(),
                             comment='Added sample_date and fraction_recovered. Fixed .15 C/K issues')
        o.metadata.change_log.append(cl)

        o.metadata.name = o.metadata.name.strip()
        if o.metadata.name[-4:].isnumeric():
            o.metadata.sample_date = o.metadata.name[-4:]
            o.metadata.name = o.metadata.name[:-5].strip()
        else:
            o.metadata.sample_date = str(o.metadata.reference.year)

        if o.oil_id == 'NO00109':  # Duplicate cuts
            o.sub_samples[0].distillation_data.cuts.pop(5)

        if o.oil_id == 'NO00120':  # non-increasing cuts
            o.sub_samples[0].distillation_data.cuts[0].vapor_temp.value = 142

        if o.oil_id in['NO00141', 'NO00141']:  # duplicate temperatures
            vis = DynamicViscosityList(o.sub_samples[2].physical_properties.dynamic_viscosities).pop(1)
            vis2 = DynamicViscosityList(o.sub_samples[2].physical_properties.dynamic_viscosities).pop(0)
            visl = DynamicViscosityList()
            visl.insert(0, vis)
            visl2 = DynamicViscosityList()
            visl2.insert(0, vis2)
            o.sub_samples[3].physical_properties.dynamic_viscosities = visl
            o.sub_samples[2].physical_properties.dynamic_viscosities = visl2

        if o.oil_id == 'NO00137':  # non-accumulate cuts
            cuts = o.sub_samples[0].distillation_data.cuts
            cuts[5].vapor_temp.value = cuts[5].vapor_temp.value + 200

        if o.oil_id == 'NO00138':  # non-accumulate cuts
            cuts = o.sub_samples[0].distillation_data.cuts
            cuts[1].fraction.value = cuts[1].fraction.value*10
            cuts[2].fraction.value = cuts[2].fraction.value*10
            cuts[3].fraction.value = cuts[3].fraction.value*10

        if o.oil_id == 'NO00130':  # Typo in Marulk oil csv
            vis = DynamicViscosityList(o.sub_samples[0].physical_properties.dynamic_viscosities).pop(1)
            visl = DynamicViscosityList()
            visl.insert(0, vis)
            o.sub_samples[0].physical_properties.dynamic_viscosities = visl

        for s in o.sub_samples:
            s.distillation_data.fraction_recovered = Concentration(max_value=1.0, unit="fraction")
            for cut in s.distillation_data.cuts:
                cut.vapor_temp.value = np.round(cut.vapor_temp.value, 1)
            pp = s.physical_properties.pour_point
            if pp is not None:
                if pp.measurement.unit == 'K':
                    off = .15
                else:
                    off = 0
                if pp.measurement.min_value is not None:
                    pp.measurement.min_value = np.floor(pp.measurement.min_value) + off
                    pp.measurement.max_value = np.floor(pp.measurement.max_value) + off
                if pp.measurement.value is not None:
                    pp.measurement.value = np.floor(pp.measurement.value) + off
            fp = s.physical_properties.flash_point
            if fp is not None:
                if fp.measurement.unit == 'K':
                    off = .15
                else:
                    off = 0
                if fp.measurement.min_value is not None:
                    fp.measurement.min_value = np.floor(fp.measurement.min_value) + off
                if fp.measurement.max_value is not None:
                    fp.measurement.max_value = np.floor(fp.measurement.max_value) + off
                if fp.measurement.value is not None:
                    fp.measurement.value = np.floor(fp.measurement.value) + off
            for dens in s.physical_properties.densities:
                if dens.ref_temp.unit == 'K':
                    off = .15
                else:
                    off = 0
                dens.ref_temp.value = np.floor(dens.ref_temp.value) + off
            for visc in s.physical_properties.dynamic_viscosities:
                if visc.ref_temp.unit == 'K':
                    off = .15
                else:
                    off = 0
                visc.ref_temp.value = np.floor(visc.ref_temp.value) + off
            for itw in s.physical_properties.interfacial_tension_water:
                if itw.ref_temp.unit == 'K':
                    off = .15
                else:
                    off = 0
                itw.ref_temp.value = np.floor(itw.ref_temp.value) + off
            for itw in s.physical_properties.interfacial_tension_seawater:
                if itw.ref_temp.unit == 'K':
                    off = .15
                else:
                    off = 0
                itw.ref_temp.value = np.floor(itw.ref_temp.value) + off

        o.status = o.validate()
        ea = len(o.validate())
        print(oil, o.metadata.name, f'({o.metadata.sample_date})', eb, ea)
        if ea > 0:
            print(o.validate())
        if folder_out is None:
            print('No writing')
        else:
            o.to_file(folder_out + o.oil_id + '.json')

if __name__ == '__main__':
    fix_norwegian_oils()
