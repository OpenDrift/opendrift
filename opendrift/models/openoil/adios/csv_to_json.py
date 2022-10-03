#! /usr/bin/env python
#
# Author: Gaute Hope (gauteh@met.no) / 2022
#
# Converts oils from old CSV format to JSON format using `adios_db` (make sure it is installed from adios_oil_database).
#
# Usage:
#
#   csv_to_json.py CSVFILE outdir/

import sys
import logging
import coloredlogs
from pathlib import Path

def main():
    from adios_db.data_sources.noaa_fm import OilLibraryCsvFile, OilLibraryRecordParser, OilLibraryAttributeMapper
    from adios_db.models.oil.validation.validate import validate_json
    from adios_db.models.oil.completeness import set_completeness

    logger = logging.getLogger(__name__)
    coloredlogs.install('debug')

    CSV = sys.argv[1]
    OUT = sys.argv[2]

    logger.info(f"Reading input file: {CSV}")
    logger.info(f"Saving JSON files to: {OUT}")

    config = None
    reader = OilLibraryCsvFile
    parser = OilLibraryRecordParser
    mapper = OilLibraryAttributeMapper

    for record in reader(CSV).get_records():
        oil_map = mapper(parser(*record))
        oil_json = oil_map.py_json()
        oil = validate_json(oil_json)
        set_completeness(oil)

        outf = Path(OUT).joinpath(f"{oil.oil_id}.json")
        logger.info(f"Writing {outf}..")
        with open(outf, 'w') as fd:
            oil.to_file(outf)

if __name__ == '__main__':
    main()

