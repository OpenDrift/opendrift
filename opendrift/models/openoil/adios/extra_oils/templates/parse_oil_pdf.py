import sys
import os
import numpy as np
import pandas as pd
import zipfile
import importlib
import adios_db.scripting as ads
from adios_db.models.oil.physical_properties import DynamicViscosityPoint

# Usage:
# Step 1:
#   parse report, and export preliminary results to a text file <oil_id>.tables:
#       $ parse_oil_pdf.py <report>.pdf >! <oil_id>.tables
#
# Step 2:
#   manually edit tables text file:
#       - delete undesired tables
#       - modify metadata, e.g. oil ID, title and authors of report
#
# Step 3:
#   use template file to create a final json file <oil_id>.json for this oil:
#       $ parse_oil_pdf.py <oil_id>.tables

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
            'add_to': 273.15,
            },
        'flash_point': {
            'names': ['Flash point', 'Flash Point'],
            'add_to': 273.15,
            },
        },
    'columns': {
        'boiling_point_temp': {'names': ['Temp.']},
        }
    }

f = sys.argv[1]  # input file name
b = f.split('.')[0]
suffix = f.split('.')[-1]
sep = ','

if suffix == 'pdf':
    # First time, the pdf is parsed and cached in a .zip file
    if not os.path.exists(b+'.zip'):
        print('Parsing ' + f)
        import camelot.io as camelot
        tables = camelot.read_pdf(f, pages='all')
        print("Total tables extracted:", tables.n)
        tables.export(b+'.csv', f='csv', compress=True)

    z = zipfile.ZipFile(b+'.zip')

    for table in z.namelist():
        df = pd.read_csv(z.open(table))
        if df.isnull().all().all():
            continue  # empty table
        if len(df.columns) not in [2, 5]:
            continue
        df = df.astype(str)
        df = df.replace(r'\n',' ', regex=True)
        headers = df.columns
        if '\n' in headers[0]:
            headers = headers[0].split('\n')
            if len(headers) > len(df.columns):
                headers = headers[0:len(df.columns)]
            elif len(headers) < len(df.columns):
                headers = df.columns
                headers = [h.replace('\n', '') for h in headers]
        maxcollengths = [df[h].astype(str).str.len().max() for h in df.columns]
        formats = ['{:>%ss}%s' % (m,sep) for m in maxcollengths]
        formats[0] = formats[0].replace('>', '<')
        formats[-1] = formats[-1].replace(sep, '')
        formats = [f.format for f in formats]
        hs = ''
        for h,f in zip(headers, formats):
            hs = hs + f(h.replace('\n', ''))
        print(hs)
        print(df.to_string(index=False, header=False, formatters=formats))
        print('\n')

    metadata = {
        'Key': 'Value',
        'oilID': 'oilID',
        'title': 'title',
        'authors': 'authors',
        'region': 'region',
        }
    
    #print('# Metadata')
    for i, (key, value) in enumerate(metadata.items()):
        print('%s%s%s' % (key, sep, value))

elif suffix == 'tables': # Step 3: parsing tables file to json
    #from pandas.compat import StringIO #if this doesn't work try: from io import StringIO
    from io import StringIO
    print('Parsing tables file')
    tables = open(f).read().split('\n\n')
    print(len(tables))
    for i, t in enumerate(tables):
        #print(t)
        #print(i)
        #title = t.split('\n')[0][2:]
        #t = '\n'.join(t.split('\n')[1:])
        df = pd.read_csv(StringIO(t), sep=sep)
        print(df, '\n')
    #print(tables, 'TAB')

stop

def parse_cuts_table(df, oil):
    pass

def parse_weathering_table(df, oil):
    max_water_cont = None
    viscosity_of_max_water_cont = None
    for i, row in df.iterrows():
        print(row, 'ROW')
        if 'weight residue' in row[0].lower():
            for c in range(1, 5):
                mo = ads.MassFraction(value=100-float(row[c]), unit="%")
                oil.sub_samples[c-1].metadata.fraction_evaporated = mo
        if 'density' in row[0].lower():
            for c in range(1, 5):
                dp = ads.DensityPoint(
                        density=ads.Density(value=float(row[c]), unit="g/ml"),
                        ref_temp=ads.Temperature(value=13, unit='C'))
                oil.sub_samples[c-1].physical_properties.densities.append(dp)
        if 'viscosity' in row[0].lower() and 'free' in row[0].lower():
            for c in range(1, 5):
                dp = DynamicViscosityPoint(
                        viscosity=ads.DynamicViscosity(value=float(row[c]), unit="cP"),
                        ref_temp=ads.Temperature(value=13, unit='C'))
                oil.sub_samples[c-1].physical_properties.dynamic_viscosities.append(dp)
        if 'viscosity' in row[0].lower() and 'max' in row[0].lower():
            viscosity_of_max_water_cont = row[1:5]
        if 'max' in row[0].lower() and 'water cont' in row[0].lower():
            max_water_cont = row[1:5]

    if viscosity_of_max_water_cont is not None and max_water_cont is not None:
        for c, visc, watcont in zip(range(1, 5), viscosity_of_max_water_cont, max_water_cont):
            try:
                eo = ads.Emulsion(water_content=ads.MassFraction(value=float(watcont), unit='%'),
                                  ref_temp=ads.Temperature(value=13, unit='C'),
                                  complex_viscosity=ads.DynamicViscosity(value=float(visc), unit='cP'))
                oil.sub_samples[c-1].environmental_behavior.emulsions.append(eo)
                print(c, eo)
            except ValueError:
                pass

    print(oil.py_json())
    oil.to_file('test.json')
    stop


#oil = ads.Oil('dummy')
#oil.sub_samples.extend([ads.Sample(), ads.Sample(), ads.Sample(), ads.Sample()])
#for i, s in enumerate(['Fresh oil', 'Topped to 150C', 'Topped to 200C', 'Topped to 250C']):
#    oil.sub_samples[i].metadata.name = s
#    oil.sub_samples[i].metadata.short_name = s


results = {}  # to store the parsed values
z = zipfile.ZipFile(b+'.zip')
print('Parsing zipfile ' + b+'.zip')
for table in z.namelist():
    #wd = {}
    #print('Parsing table ' + table)
    df = pd.read_csv(z.open(table))
    if df.isnull().all().all():
        continue  # empty table
    print('# ' + table)
    print('%s\n' % df.to_string())
    continue
    if table != 'TOR_oil-page-95-table-1.csv':
        continue

    parse_weathering_table(df, oil)
    parsing

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
    stop

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
