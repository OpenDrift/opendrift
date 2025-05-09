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

def open_dataset_opendrift(source, zarr_storage_options=None, open_mfdataset_options={}, chunks=None):
    """ Wrapper around Xarray open_dataset and open_mfdataset.

    Common wrapper/opener to be used for all Xarray based readers

    xarray.open_dataset will be used if source is:
        - a single netCDF file or OPeNDAP URL

    xarray.open_mfdataset will be used if source is:
        - a list of netCDF files
        - a filename with wildcards (* ? or [)

    cf-times are decoded after removing any offending variables
        (e.g. if units equals "hours since analysis")

    """

    if chunks is not None:
        open_mfdataset_options['chunks'] = chunks

    if isinstance(source, xr.Dataset):
        ds = source
    elif zarr_storage_options is not None:  # This could better be handled outside of this method
        ds = xr.open_zarr(source, storage_options=zarr_storage_options)
        ds.name = source
    elif isinstance(source, list) or any(s in str(source) for s in [ '*', '?', '[' ]):
        logger.info('Opening files with xarray.open_mfdataset')
        ds = xr.open_mfdataset(source, data_vars='minimal', coords='minimal', compat='override',
                               decode_times=False, **open_mfdataset_options)
    else:
        logger.info('Opening file with xr.open_dataset')
        ds = xr.open_dataset(source, decode_times=False)

    # Decode CF times
    offending = ds.filter_by_attrs(units='hours since analysis')  # Found e.g. in HYCOM datasets
    if len(offending) > 0:
        logger.warning(f'Removing variables that cannot be CF decoded: {list(offending.variables)}')
        ds = ds.drop_vars(offending)
    ds = xr.decode_cf(ds, decode_times=True)

    # Chunk size of time dimension should be 1
    ds = ds.unify_chunks()
    for dim,chunksize in ds.chunks.items():
        chunksize = chunksize[0]
        if 'time' in dim:
            if chunksize > 1:
                # Note: ds = ds.chunk({dim: 1}) does not have the desired effect
                logger.warning(f'Chunk size for dimension {dim} is {chunksize}. Reopening dataset with chunk {dim}:1.')
                ds = open_dataset_opendrift(source, zarr_storage_options=zarr_storage_options,
                                            open_mfdataset_options=open_mfdataset_options, chunks={dim: 1})

    return ds

def datetime_from_variable(var):
    import pandas as pd
    try:
        return pd.to_datetime(var).to_pydatetime()
    except:
        logger.warning('Could not decode time with Pandas')
        datetimeindex = var.to_index().to_datetimeindex()
        times = pd.to_datetime(datetimeindex).to_pydatetime()
        logger.info('Decoded time through datetimeindex')
        return times

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

def applicable_readers(url):
    '''Return a list of readers that are possible candidates for a given URL, filename or product ID'''

    from opendrift.readers import reader_netCDF_CF_generic
    from opendrift.readers import reader_netCDF_CF_unstructured
    from opendrift.readers import reader_ROMS_native
    from opendrift.readers import reader_copernicusmarine
    if len(glob.glob(url)) > 0 or any(e in url for e in [':', '/']):
        return [reader_netCDF_CF_generic, reader_ROMS_native, reader_netCDF_CF_unstructured]
    elif '_' in url or '-' in url:  # should have better indentificator
        return [reader_copernicusmarine] 
    else:
        return []

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

    reader_modules = applicable_readers(url)

    for reader_module in reader_modules:
        try:
            logger.debug(f'Testing reader {reader_module}')
            r = reader_module.Reader(url)
            return r
        except Exception as e:
            logger.debug('Could not open %s with %s' % (url, reader_module))

    return None  # No readers worked
