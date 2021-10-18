#!/usr/bin/env python

from opendrift.models.leeway import Leeway
from opendrift.models.openoil import OpenOil
from opendrift.models.shipdrift import ShipDrift
from datetime import datetime, timedelta
import click
import sys
import json
import math


def get_leeway_model(loglevel):
    model = Leeway(loglevel=loglevel)
    set_config_override(model)
    return model

def get_openoil_model(loglevel):
    model = OpenOil(location='NORWAY', loglevel=loglevel)
    set_config_override(model)
    return model

def get_shipdrift_model(loglevel):
    model = ShipDrift(loglevel=loglevel)
    set_config_override(model)
    return model


def set_config_override(model):
    model._config['seed:number']['default'] = 1000
    model._config['seed:number']['value'] = 1000
    model._config['seed:number']['max'] = 100000

    model._config['general:time_step_minutes']['default'] = 15
    model._config['general:time_step_minutes']['value'] = 15
    model._config['general:time_step_minutes']['min'] = 5

    model._config['general:time_step_output_minutes']['default'] = 60
    model._config['general:time_step_output_minutes']['value'] = 60
    model._config['general:time_step_output_minutes']['min'] = 30

    try:
        model._config['seed:seafloor']['level'] = 2
        model._config['seed:z']['level'] = 2
    except KeyError:
        pass

models = {
    'leeway': get_leeway_model,
    'openoil': get_openoil_model,
    'shipdrift': get_shipdrift_model
}


plots = {
    'plot': lambda model, filename: model.plot(filename=filename),
    'heatmap': lambda model, filename: model.write_netcdf_density_map(filename=filename),
    'text': lambda model, filename: model.export_ascii(filename),
    'oil_budget': lambda model, filename: model.plot_oil_budget(filename=filename)
}


@click.command()
@click.option('-m', '--model', required=True, type=click.Choice(list(models.keys())), help="run the given model")
@click.option('-o', '--output', type=click.Path(exists=False, writable=True), help='Plot to the given file')
@click.option('--latitude', type=float, multiple=True, help='Use the given latitudes')
@click.option('--longitude', type=float, multiple=True, help='Use the given longitudes')
@click.option('--radius', type=int, multiple=True, help='Seed with the given radius')
@click.option('--time', multiple=True, help='Seed at the given time')
@click.option('--polygon', help='Seed with the given polygon instead of latitude/longitude. Format like this: "4.7 60.34, 4.7 60.4, 4.8 60.4, 4.8 60.34, 4.7 60.34"')
@click.option('-c', '--config', help='Set the given config option', multiple=True, type=str)
@click.option('-d', '--forcing_data', multiple=True, help='Use the given reader. This option can be used several times to use several readers.')
@click.option('--duration', default=6, type=int, help='run simulation for th given number of hours. Negative values are allowed, which will cause the simulation to run backwards')
@click.option('--plot', multiple=True, help='Create the additional plots. Use syntax "--plot heatmap=file.nc". Available plots: ' + ','.join(plots.keys()) + '.')
@click.option('--list-config', default=False, help='List available config options instead of running a simulation', is_flag=True)
@click.option('--print-summary', default=False, help='End simulation by printing a summary of timings and forcing data that wass used', is_flag=True)
@click.option('--loglevel', type=int, default=40, help='Logging verbosity 0 is everything, 50 is nothing')
def main(model, output, latitude, longitude, radius, time, polygon, config, forcing_data, duration, plot, list_config, print_summary, loglevel):
    '''Run an opendrift simulation. Use latitude/longitude/radius/time as many times as you wish to spread seeds out over a line.'''
    if model not in models:
        raise Exception('model must be specified')

    if list_config:
        display_config_options(models[model])
        sys.exit(0)

    if not output:
        raise Exception('missing output file')

    if not forcing_data:
        raise Exception('no forcing data provided')

    if not polygon:
        if not latitude:
            raise Exception(
                'missing location spec (latitude, longitude, radius and time)')
        if len(longitude) != len(latitude) or len(latitude) != len(radius) or len(radius) != len(time):
            raise Exception(
                'latitude, longitude, radius and time must have the same number of elements')
    elif latitude or longitude or radius:
        raise Exception('cannot mix latitude/longitude/radius with polygon')
    elif len(time) != 1:
        raise Exception('when using polygon, you must have exactly one time parameter')


    model = models[model](loglevel=loglevel)

    model.add_readers_from_list(forcing_data)

    available_config = model.get_configspec()
    for cfg in config:
        key, val = cfg.split('=', 1)
        try:
            config_spec = available_config[key]
            conf_type = config_spec['type']
            if conf_type == 'int':
                val = int(val)
            elif conf_type == 'float':
                val = float(val)
            elif conf_type == 'bool':
                val = val.lower() == 'true'
        except:
            raise Exception('invalid configuration key: ' + key)
        model.set_config(key, val)

    if polygon:
        model.seed_from_wkt(
            wkt = 'MULTIPOLYGON(((' + polygon + ')))',
            time=datetime.strptime(time[0], "%Y-%m-%dT%H:%M:%SZ")
        )
    else:
        if len(longitude) == 1:
            model.seed_elements(
                lon=longitude[0],
                lat=latitude[0],
                radius=radius[0],
                time=datetime.strptime(time[0], "%Y-%m-%dT%H:%M:%SZ")
            )
        else:
            model.seed_cone(
                lon=longitude,
                lat=latitude,
                radius=radius,
                time=[datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ") for t in time],
                number=1000
            )

    if duration < 0:
        # If duration is negative, time_step_minutes and time_step_output_minutes must alos be negative
        time_step_minutes = model._config['general:time_step_minutes']['value']
        if time_step_minutes > 0:
            model._config['general:time_step_minutes']['value'] = -time_step_minutes

        time_step_output_minutes = model._config['general:time_step_output_minutes']['value']
        if time_step_output_minutes and time_step_output_minutes > 0:
            model._config['general:time_step_output_minutes']['value'] = -time_step_output_minutes

    model.run(
        outfile=output, 
        duration=timedelta(hours=duration)
    )

    for plot in plot:
        plotType, filename = plot.split('=', 2)
        plots[plotType](model, filename)

    if print_summary:
        print(model, file=sys.stderr)


def display_config_options(model):
    model = model(loglevel=50)
    spec = model.get_configspec()
    for name, s in spec.items():
        if 'max' in s and s['max'] == math.inf:
            del(s['max'])
        if 'min' in s and s['min'] == math.inf:
            del(s['min'])
    doc = json.dumps(spec, indent=2, allow_nan=False)
    click.echo(doc)


if __name__ == '__main__':
    main()
