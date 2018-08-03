import logging
import glob
from opendrift.readers.reader_netCDF_CF_generic import Reader

def reader_from_url(url, timeout=10):
    '''Make readers from URLs or paths to datasets'''

    if isinstance(url, list):
        return [reader_from_url(u) for u in url]

    files = glob.glob(url)
    for f in files:  # Regular file
        try:
            r = Reader(f)
            return r

        except:
            logging.warning('%s is not a netCDF CF file recognised by '
                            'OpenDrift' % f)
            try:
                from opendrift.readers.reader_ROMS_native import Reader as Reader_ROMS_native
                r = Reader_ROMS_native(f)
                return r
            except:
                logging.warning('%s is also not a ROMS netCDF file recognised by '
                                'OpenDrift' % f)
                try:
                    from opendrift.readers.reader_grib import Reader as Reader_grib
                    r = Reader_grib(f)
                    return r
                except:
                    logging.warning('%s is also not a GRIB file recognised by '
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
                logging.debug('Applied NETRC credentials')
            except:
                logging.debug('Could not apply NETRC credentials')
            urllib_request.urlopen(request, timeout=timeout)
        except Exception as e:
            # Error code 400 is expected!
            if not isinstance(e, urllib_request.HTTPError) or e.code != 400:
                logging.warning('ULR %s not accessible: ' % url + str(e))
                return None
            try:
                r = Reader(url)
                return r
            except Exception as e:
                logging.warning('%s is not a netCDF file recognised '
                                'by OpenDrift: %s' % (url, str(e)))
                return None
