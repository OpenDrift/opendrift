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

import importlib
import logging; logger = logging.getLogger(__name__)
import glob
import json

def reader_from_url(url, timeout=10):
    '''Make readers from URLs or paths to datasets'''

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

    from opendrift.readers import reader_netCDF_CF_generic

    if isinstance(url, list):
        return [reader_from_url(u) for u in url]

    files = glob.glob(url)
    for f in files:  # Regular file
        try:
            r = reader_netCDF_CF_generic.Reader(f)
            return r

        except:
            logger.warning('%s is not a netCDF CF file recognised by '
                            'OpenDrift' % f)
            try:
                from opendrift.readers.reader_ROMS_native import Reader as Reader_ROMS_native
                r = Reader_ROMS_native(f)
                return r
            except:
                logger.warning('%s is also not a ROMS netCDF file recognised by '
                                'OpenDrift' % f)
                try:
                    from opendrift.readers.reader_grib import Reader as Reader_grib
                    r = Reader_grib(f)
                    return r
                except:
                    logger.warning('%s is also not a GRIB file recognised by '
                                    'OpenDrift' % f)

    if files == []:  # Try with OPeNDAP URL
        try:  # Check URL accessibility/timeout
            try:
                # for python 3
                import urllib.request as urllib_request
            except ImportError:
                # for python 2
                import urllib2 as urllib_request
            request = urllib_request.Request(url)
            try:  # netrc
                import netrc
                import base64
                parts = urllib_request.urlparse.urlparse(url)
                login, account, password = netrc.netrc().authenticators(parts.netloc)
                creds = base64.encodestring('%s:%s' % (login, password)).strip()
                request.add_header("Authorization", "Basic %s" % creds)
                logger.debug('Applied NETRC credentials')
            except:
                logger.debug('Could not apply NETRC credentials')
            urllib_request.urlopen(request, timeout=timeout)
        except Exception as e:
            # Error code 400 is expected!
            if not isinstance(e, urllib_request.HTTPError) or e.code != 400:
                status_code = None
                try:  # Trying with requests library (should be default)
                    import requests
                    resp = requests.get(url + '.das')
                    status_code = resp.status_code
                except Exception as e:
                    logger.warning('ULR %s not accessible: ' % url + str(e))
                    return None

                if status_code >= 400:
                    logger.warning('ULR %s not accessible: ' % url + str(e))
                    return None
            try:
                from opendrift.readers import reader_netCDF_CF_generic
                r = reader_netCDF_CF_generic.Reader(url)
                return r
            except Exception as e:
                logger.warning('%s is not a netCDF file recognised '
                                'by OpenDrift: %s' % (url, str(e)))
                try:
                    from opendrift.readers.reader_ROMS_native import Reader as Reader_ROMS_native
                    r = Reader_ROMS_native(url)
                    return r
                except Exception as e:
                    logger.warning('%s is also not a ROMS netCDF file recognised by '
                                    'OpenDrift: %s' % (url, str(e)))

                return None
