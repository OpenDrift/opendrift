import os
import numpy as np
import pandas as pd
import zipfile

template = {
    'rows': {
        'density': {
            'names': ['Specific gravity'],
            'factor': 1000  # Convert g/l to kg/m3
            },
        'weathered': {
            'names': ['Residue (wt'],
            'subtract_from': 100
            },
        #'viscosity': {
        #    'names': ['Viscosity at 13'],
        #    'factor': 1
        #    },
        'wax_content': {
            'names': ['Wax Content'],
            'factor': 0.01
            },
        'asphaltenes': {
            'names': ['Asphaltenes'],
            'factor': 0.01
            },
        'pour_point': {
            'names': ['Pour Point'],
            },
        },
    'columns': {
        'boiling_point_temp': {'names': ['Temp.(']},
        }
    }

f = 'EDOCS-#51667-v1-Njord_2003_Forvitringsrapport.PDF'
b = f.split('.')[0]

if not os.path.exists(b+'.zip'):
    print('Parsing ' + f)
    import camelot
    tables = camelot.read_pdf(f, pages='all')
    print("Total tables extracted:", tables.n)
    tables.export(b+'.csv', f='csv', compress=True)


results = {}  # to store the parsed values
z = zipfile.ZipFile(b+'.zip')
print('Parsing zipfile ' + b+'.zip')
for table in z.namelist():
    #print('Parsing table ' + table)
    df = pd.read_csv(z.open(table))
    # Search for row parameters
    for v, r in template['rows'].items():
        for name in r['names']:
            for i, row in df.iterrows():
                if name in str(row):
                    results[v] = []
                    try:
                        for aval in np.array([rva.replace(',', '.') for rva in np.array(row.values)]):
                            try:
                                results[v].append(np.float32(aval.replace(',', '.')))
                            except:
                                pass
                        results[v] = np.array(results[v])
                        if 'factor' in r:
                            results[v] = results[v]*r['factor']
                        if 'subtract_from' in r:
                            results[v] = r['subtract_from'] - results[v]
                    except:
                        #print('Could not parse row: ' + str(row))
                        pass

    for v, r in template['columns'].items():
        for name in r['names']:
            for i, row in df.iterrows():
                if name in str(row):
                    print(df)
                    #jakka

for n, r in results.items():
    print(n, r)
