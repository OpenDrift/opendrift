To add new oils:

- run parser: openoil/adios/extra_oils/templates/parse_oil_pdf.py <PDF-report>
- copy and edit new py-file in openoil/adios/extra_oils/python_oils
- run parse.py in the same folder to create json file
- copy json file to software/noaa-oil-data/data/oil/NO
- Commit and push to branch production
