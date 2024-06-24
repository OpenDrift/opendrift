"""
Readers
=======

Readers are responsible for providing Opendrift with data about the
`enviornment` of the drifting particles or elements.

All readers descend from :class:`.basereader.BaseReader`. A reader generally also descends from one of the few general classes of readers. When writing a new reader consider which one fits your input data best:

    * :class:`.basereader.continuous.ContinuousReader`
    * :class:`.basereader.structured.StructuredReader`
    * :class:`.basereader.unstructured.UnstructuredReader`

The `ContinuousReader` is suited for data that can be defined at any point within the domain, or if the reader does its own interpolation internally. E.g. a :class:`synthetic eddy <.reader_ArtificialOceanEddy.Reader>`, or a :class:`constant <.reader_constant.Reader>`. The `StructuredReader` aids in interpolation when creating a reader of data on a :class:`regular grid <.reader_netCDF_CF_generic.Reader>`, while the `UnstructuredReader` provides the basics for :class:`irregularily gridded data <.reader_netCDF_CF_unstructured.Reader>` (e.g. finite volume models).

.. seealso::

    See the reader-types or reader-implementations for more details.

    See :class:`.basereader.BaseReader` for how readers work internally.
"""

from datetime import datetime, timedelta
import importlib
import logging; logger = logging.getLogger(__name__)
import glob
import json
import opendrift
import xarray as xr

def open_mfdataset_overlap(url_base, time_series=None, start_time=None, end_time=None, freq=None, timedim='time'):
    if time_series is None:
        construct_from_times
    urls = [t.strftime(url_base) for t in time_series]
    time_step = time_series[1] - time_series[0]
    print('Opening individual URLs...')
    chunks = {'time': 1, 'depth': 1, 'Y': 17, 'X': 2602}
    datasets = [xr.open_dataset(u, chunks=chunks).sel({timedim: slice(t, t+time_step-timedelta(seconds=1))})
                for u,t in zip(urls, time_series)]
    print('Concatenating...')
    ds = xr.concat(datasets, dim=timedim,
                   compat='override', combine_attrs='override', join='override', coords='minimal', data_vars='minimal')
    return ds

def reader_from_url(url, timeout=10):
    '''Make readers from URLs or paths to datasets'''

    if isinstance(url, list):
        return [reader_from_url(u) for u in url]

    try:  # Initialise reader from JSON string
        j = json.loads(url)
        try:
            reader_module = importlib.import_module(
                    'opendrift.readers.' + j['reader'])
            reader = getattr(reader_module, 'Reader')
            del j['reader']
            reader = reader(**j)
            return reader
        except Exception as e:
            logger.warning('Creating reader from JSON failed:')
            logger.warning(e)
    except:
        pass

    if len(glob.glob(url)) == 0:  # Check if this is a URL, and not giving timeout
        import requests
        try:
            resp = requests.get(url, timeout=timeout)
        except requests.exceptions.MissingSchema:

            logger.info('Neither a file or a URL, could still be a CMEMS product ID: ' + url)
            #return None
        except:
            logger.info('Connection error for ' + url)
            return None

    reader_modules = ['reader_netCDF_CF_generic',
                      'reader_ROMS_native',
                      'reader_copernicusmarine']

    for rm in reader_modules:
        reader_module = importlib.import_module('opendrift.readers.' + rm)
        try:
            r = reader_module.Reader(url)
            return r
        except Exception as e:
            print('Could not open %s with %s' % (url, rm))

    return None  # No readers worked
