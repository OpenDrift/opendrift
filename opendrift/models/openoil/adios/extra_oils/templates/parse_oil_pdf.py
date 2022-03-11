import sys
import os
import numpy as np
import pandas as pd
import zipfile

np.set_printoptions(suppress=True)

template = {
    'rows': {
        'density': {
            'names': ['Specific gravity', 'Density'],
            'factor': 1000  # Convert g/l to kg/m3
            },
        'weathered': {
            'names': ['Residue (wt', 'Weight Residue', 'Vol, Topped'],
            #'names': ['Weight Residue'],
            'subtract_from': 1,
            'factor': .01
            },
        'viscosity': {
            'names': ['Viscosity at 13°C (cP) to predict', 'Viscosity of water-free'],
            'factor': 0.001
            },
        'wax_content': {
            'names': ['Wax Content'],
            'factor': 0.01
            },
        'asphaltenes': {
            'names': ['Asphaltenes'],
            'factor': 0.01
            },
        'pour_point': {
            'names': ['Pour point', 'Pour Point'],
            'add_to': 273,
            },
        'flash_point': {
            'names': ['Flash point', 'Flash Point'],
            'add_to': 273,
            },
        },
    'columns': {
        'boiling_point_temp': {'names': ['Temp.']},
        }
    }

f = sys.argv[1]
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
    #print(df)
    # Search for row parameters
    for v, r in template['rows'].items():
        for name in r['names']:
            for i, row in df.iterrows():
                if name in str(row):
                    varname = v
                    while(varname) in results:
                        varname = varname + 'a'
                    results[varname] = []
                    #print(varname + str(row.values))
                    for ri in range(len(row.values)):
                        try:
                            row.values[ri] = np.float32(row.values[ri].replace(',', '.'))
                        except:
                            pass
                    try:
                        for aval in row.values:
                            try:
                                results[varname].append(np.float32(aval))
                            except:
                                pass
                        results[varname] = np.array(results[varname])
                        if 'factor' in r:
                            results[varname] = results[varname]*r['factor']
                        if 'subtract_from' in r:
                            results[varname] = r['subtract_from'] - results[varname]
                        if 'add_to' in r:
                            results[varname] = r['add_to'] + results[varname]
                        #print(results[varname])
                    except:
                        print('Could not parse row: ' + str(row))
                        stop
                        pass

    #print('='*30)
    #print(df)
    #print('='*30)
    for v, r in template['columns'].items():
        headerline = df.columns
        if 'Temp' in headerline[0] and ('Vol' in headerline[1] or 'Duva' in headerline[1]):
            cut_temps = df.iloc[:,0].values + 273
            cut_vols = df.iloc[:,1].values
            try:
                cut_vols = cut_vols / 100
            except:
                cut_vols = [np.float64(val.replace(',', '.'))/100 for val in cut_vols]

            re = np.empty(len(cut_temps)*3)
            for n in range(len(cut_temps)):
                re[n*3] = cut_temps[n]
                re[n*3+2] = cut_vols[n]
            print('Cuts:\n' + np.array2string(re.astype(np.float32), separator='\t').replace('\n', '\t'))
#print('='*30)
#print(pd.read_csv(z.open(z.namelist()[0])))
#print('='*30)
for n, r in results.items():
    if len(r) > 0 and np.sum(np.isfinite(r))>0:
        print(n, r)
print('='*30)

for aas in ['', 'a', 'aa', 'aaa', 'aaaa', 'aaaaa']:
    dstring = 'density' + aas
    vstring = 'viscosity' + aas
    wstring = 'weathered' + aas
    rd = np.empty(12)
    rv = np.empty(12)
    for i in range(4):
        if dstring in results and len(results[dstring])>0:
            rd[i*3] = results[dstring][i]
            rd[i*3+1] = 286
            rd[i*3+2] = results[wstring][i]
        if vstring in results and len(results[vstring])>0:
            rv[i*3] = results[vstring][i]
            rv[i*3+1] = 286
            rv[i*3+2] = results[wstring][i]

    if vstring in results:
        print(vstring, np.array2string(rv.astype(np.float32), separator='\t').replace('\n', '\t'))
    if dstring in results:
        print(dstring, np.array2string(rd.astype(np.float32), separator='\t').replace('\n', '\t'))

    
#try:
#    r = np.empty(12)
#    for i in range(4):
#        r[i*3] = results['density'][i]
#        r[i*3+1] = 286
#        r[i*3+2] = results['weathered'][i]
#    print(np.array2string(r.astype(np.float32), separator='\t').replace('\n', '\t'))
#except:
#    print('No density')
#
#try:
#    for i in range(4):
#        r[i*3] = results['viscosity'][i]
#    print(np.array2string(r.astype(np.float32), separator='\t').replace('\n', '\t'))
#except:
#    print('No viscosity')
