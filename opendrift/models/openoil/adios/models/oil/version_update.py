"""
version_update.py

code for dealing with different versions of the data model

NOTE: this may need some refactoring when it gets more complicated

"""

from .version import Version, VersionError

from .oil import ADIOS_DATA_MODEL_VERSION


class Updater:
    ver_from = None
    ver_to = None

    def __call__(self, py_json):
        this_ver = Version(py_json.get('adios_data_model_version'))
        if this_ver == self.ver_to:  # nothing to be done
            return py_json
        elif this_ver != self.ver_from:  # can't update this
            return py_json
        else:  # we can do the update
            return self.update(py_json)

class update_0_10_0_to_0_10_0(Updater):

    ver_from = Version(0, 10, 0)
    ver_to = Version(0, 11, 0)

    def update(self, py_json):
        # just in case
        this_ver = Version(py_json.get('adios_data_model_version'))

        if this_ver != self.ver_from:
            raise ValueError("Update called with JSON of wrong version"
                             f"JSON version: {this_ver}, this updater "
                             f"can update: {self.ver_from}")
        # change the name of the fraction_weathered attribute
        if 'sub_samples' in py_json:  # very sparse records may not
            for ss in py_json['sub_samples']:
                md = ss['metadata']
                # note: this may add a fraction_weathered = None, but that's OK
                md['fraction_evaporated'] = md.get('fraction_weathered')
                md.pop('fraction_weathered', None)
        py_json['adios_data_model_version'] = str(self.ver_to)

        return py_json


# NOTE: updaters need to be in order
#       not hard when there is only one :-)
UPDATERS = [update_0_10_0_to_0_10_0()]


def update_json(py_json):
    """
    updates JSON for an oil object from an older version to a newer one
    """

    cur_ver = ADIOS_DATA_MODEL_VERSION
    try:
        ver = Version(py_json['adios_data_model_version'])
    except KeyError:
        # assume it's the version from before we added a version
        ver = Version(0, 10, 0)
        py_json['adios_data_model_version'] = str(ver)
    if ver == cur_ver:  # nothing to be done
        return py_json
    # run it through the updaters
    for update in UPDATERS:
        py_json = update(py_json)
    # Now see if it worked
    ver = Version(py_json.get('adios_data_model_version'))
    if ver == cur_ver:
        return py_json
    elif ver > cur_ver:
        raise VersionError(f"Version: {ver} is not supported by this version of Oil object")
    else:
        raise VersionError(f"updater not available for version: {ver}")
