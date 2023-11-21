import sys
import os
import numpy as np
import pandas as pd
import zipfile
import importlib
import traceback
import tkinter as tk
from tkinter import ttk

def parse_weathering_table(df, oil):
    import adios_db.scripting as ads
    from adios_db.models.oil.physical_properties import DynamicViscosityPoint

    max_water_cont = None
    viscosity_of_max_water_cont = None
    try:
        for i, row in df.iterrows():
            print(row, 'ROW')
            if 'weight residue' in row[0].lower():
                for c in range(1, 5):
                    mo = ads.MassFraction(value=100-float(row[c]), unit="%")
                    oil.sub_samples[c-1].metadata.fraction_evaporated = mo
            elif 'density' in row[0].lower():
                for c in range(1, 5):
                    dp = ads.DensityPoint(
                            density=ads.Density(value=float(row[c]), unit="g/ml"),
                            ref_temp=ads.Temperature(value=13, unit='C'))
                    oil.sub_samples[c-1].physical_properties.densities.append(dp)
            elif 'viscosity' in row[0].lower() and 'free' in row[0].lower():
                for c in range(1, 5):
                    dp = DynamicViscosityPoint(
                            viscosity=ads.DynamicViscosity(value=float(row[c]), unit="cP"),
                            ref_temp=ads.Temperature(value=13, unit='C'))
                    oil.sub_samples[c-1].physical_properties.dynamic_viscosities.append(dp)
            elif 'viscosity' in row[0].lower() and 'max' in row[0].lower():
                viscosity_of_max_water_cont = row[1:5]
            elif 'max' in row[0].lower() and 'water cont' in row[0].lower():
                max_water_cont = row[1:5]
    except:
        raise ValueError('Not a weathering table')

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
    else:
        logfile.write('Discarding whole table\n')
        return

    print(oil, 'Finished oil!')

def parse_cuts_table(df, oil):
    print('PARSING CUTS')
    print(df)
    cd = df.to_numpy()
    print(cd, cd.shape)
    temperatures = cd[:,0]
    volumes = cd[:,1]
    print('PARSED CUTS')

class ScrollableFrame(tk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = tk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        container.pack()
        canvas.pack(side="left", fill="x", expand=0)
        scrollbar.pack(side="right", fill="y")

class ParseOilPDF(tk.Tk):

    def __init__(self):
        f = sys.argv[1]  # input file name
        b = f.split('.')[0]

        if not os.path.exists(b+'.zip'):
            print('Parsing %s, this may take some minutes...' % f)
            import camelot.io as camelot
            tables = camelot.read_pdf(f, pages='all', backend='poppler')
            print('Total tables extracted: %s\n' % tables.n)
            tables.export(b+'.csv', f='csv', compress=True)

        z = zipfile.ZipFile(b+'.zip')
        print('Number of tables: %s\n' % len(z.namelist()))

        tk.Tk.__init__(self)
        self.title('Oil parser for Sintef PDF reports')

        # Make tabs
        self.n = ttk.Notebook(self.master)
        self.n.pack()
        self.tabs = ttk.Frame(self.n)
        self.oil_main = ttk.Frame(self.n)
        self.oil = ScrollableFrame(self.oil_main)
        self.n.add(self.tabs, text='Tables')
        self.n.add(self.oil_main, text='Oil')

        mf = tk.Frame(self.tabs, bg='cyan', width=1500, height=150, padx=3, pady=3)
        cf = tk.Frame(self.tabs, bg='yellow', width=1500, height=1500, padx=3, pady=3)

        mf.pack()
        cf.pack(fill=tk.BOTH, expand=1)

        properties = ScrollableFrame(cf)
        cuts = ScrollableFrame(cf)
        weathering = ScrollableFrame(cf)

        tk.Label(properties.scrollable_frame, bg='red', text="Oil properties",
                 font=('Courier', 18, 'bold')).pack()
        tk.Label(cuts.scrollable_frame, bg='red', text="Cuts",
                  font=('Courier', 18, 'bold')).pack()
        tk.Label(weathering.scrollable_frame, bg='red', text="Weathering tables",
                  font=('Courier', 18, 'bold')).pack()

        self.tables = {}
        self.dataframes = {}
        self.selected_tables = []
        tableno = 0
        for table in z.namelist():
            df = pd.read_csv(z.open(table))
            if df.iloc[:, 1:].isnull().all().all():
                print('Skipping empty table\n')
                continue  # Do not print empty tables
            if len(df.columns) not in [2, 5]:
                print('Skipping table without 2 or 5 columns\n')
                continue

            color='lightgray'
            if len(df.columns) == 2:
                if 'temp' in df.columns[0].lower() and (
                        'vol' in df.columns[0].lower() or 'vol' in df.columns[1].lower()):
                    parent = cuts
                    width=60
                else:
                    parent = properties
                    width=60
            elif len(df.columns) == 5:
                # Check if this is a valid weathering table
                oil = ads.Oil('dummy')
                oil.sub_samples.extend([ads.Sample(), ads.Sample(), ads.Sample(), ads.Sample()])
                for i, s in enumerate(['Fresh oil', 'Topped to 150C', 'Topped to 200C', 'Topped to 250C']):
                    oil.sub_samples[i].metadata.name = s
                    oil.sub_samples[i].metadata.short_name = s
                try:
                    parse_weathering_table(df, oil)
                except Exception as e:
                    print(e)
                    print(traceback.format_exc())
                    print('Not a weathering table')
                    continue

                parent = weathering
                width=150

            self.dataframes[table] = df
            txt = str(df).split('\n')
            width = len(txt[0])
            for tx in txt:
                print(tx, len(tx), 'JOLA')
            color='yellow'
            self.tables[table] = tk.Text(parent.scrollable_frame, bg=color, width=width)
            tableno = tableno + 1
            self.tables[table].insert(tk.END, df)
            self.tables[table].pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)
            self.tables[table].bind('<Button-1>', lambda event,s=table:self.clicktable(s))

        properties.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        cuts.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        weathering.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

    def clicktable(self, table):
        if table in self.selected_tables:  # Remove/unselect
            self.tables[table].config(bg='lightgray')
            self.selected_tables.remove(table)
        else:
            self.tables[table].config(bg='green')
            self.selected_tables.append(table)

        print(self.selected_tables)
        print('Making oil')
        oil = ads.Oil('dummy')
        oil.sub_samples.extend([ads.Sample(), ads.Sample(), ads.Sample(), ads.Sample()])
        for i, s in enumerate(['Fresh oil', 'Topped to 150C', 'Topped to 200C', 'Topped to 250C']):
            oil.sub_samples[i].metadata.name = s
            oil.sub_samples[i].metadata.short_name = s
        print(oil)
        for tablename in self.selected_tables:
            df = self.dataframes[tablename]
            print('-'*22)
            print(df, 'DF')
            print('-'*22)
            try:
                parse_weathering_table(df, oil)
            except Exception as e:
                print(e)
                print(traceback.format_exc())
            try:
                parse_cuts_table(df, oil)
            except Exception as e:
                print(e)
                print(traceback.format_exc())

        oil.to_file('test.json')
        v = oil.validate()
        print(v, 'VALIDATE')
        import os
        os.system('cat test.json')
        with open('test.json', 'r') as file:
            oil_json = file.read()#.replace('\n', '')
        for widget in self.oil.scrollable_frame.winfo_children():  # Clear oil tab
            widget.destroy()
        tk.Label(self.oil, text=oil_json,
                 font=('Courier', 10, 'normal')).pack()

if __name__ == '__main__':

    import adios_db.scripting as ads

    ParseOilPDF().mainloop()


########################################################
## Step 1: the pdf is parsed and cached in a .zip file
########################################################
#if suffix == 'pdf':
#    if not os.path.exists(b+'.zip'):
#        logfile.write('Parsing ' + f)
#        import camelot.io as camelot
#        tables = camelot.read_pdf(f, pages='all')
#        logfile.write('Total tables extracted: %s\n' % tables.n)
#        tables.export(b+'.csv', f='csv', compress=True)
#
#    z = zipfile.ZipFile(b+'.zip')
#
#    for table in z.namelist():
#        df = pd.read_csv(z.open(table))
#        logfile.write('*'*33 + '\n')
#        logfile.write(str(df) + '\n')
#        logfile.write('-'*33 + '\n')
#        if df.iloc[:, 1:].isnull().all().all():
#            logfile.write('Skipping empty table\n')
#            continue  # Do not print empty tables
#        if len(df.columns) not in [2, 5]:
#            logfile.write('Skipping table without 2 or 5 columns\n')
#            continue
#        df = df.astype(str)
#        df = df.replace(r'\n',' ', regex=True)
#        headers = df.columns
#        if '\n' in headers[0]:
#            headers = headers[0].split('\n')
#            if len(headers) > len(df.columns):
#                headers = headers[0:len(df.columns)]
#            elif len(headers) < len(df.columns):
#                headers = df.columns
#                headers = [h.replace('\n', '') for h in headers]
#        maxcollengths = [df[h].astype(str).str.len().max() for h in df.columns]
#        formats = ['{:>%ss}%s' % (m,sep) for m in maxcollengths]
#        formats[0] = formats[0].replace('>', '<')
#        formats[-1] = formats[-1].replace(sep, '')
#        formats = [f.format for f in formats]
#        hs = ''
#        for h,f in zip(headers, formats):
#            hs = hs + f(h.replace('\n', ''))
#        logfile.write(str(hs) + '\n')
#        outfile.write(df.to_string(index=False, header=True, formatters=formats) + '\n\n')
#        logfile.write('*'*33 + '\n')
#
#    metadata = {
#        'Key': 'Value',
#        'oilID': 'oilID',
#        'title': 'title',
#        'authors': 'authors',
#        'region': 'region',
#        }
#
#    #print('# Metadata')
#    for i, (key, value) in enumerate(metadata.items()):
#        outfile.write('%s%s%s\n' % (key, sep, value))
#    outfile.write('\n')
#
###########################################
## Step 3: parsing tables file to create
##         oil object, and export to json
###########################################
#elif suffix == 'tables':
#
#    logfile.write('='*33 + '\n')
#    oil = ads.Oil('dummy')
#    oil.sub_samples.extend([ads.Sample(), ads.Sample(), ads.Sample(), ads.Sample()])
#    for i, s in enumerate(['Fresh oil', 'Topped to 150C', 'Topped to 200C', 'Topped to 250C']):
#        oil.sub_samples[i].metadata.name = s
#        oil.sub_samples[i].metadata.short_name = s
#
#    #from pandas.compat import StringIO #if this doesn't work try: from io import StringIO
#    from io import StringIO
#    logfile.write('Parsing tables file\n')
#    tables = open(f).read().split('\n\n')
#    logfile.write('Number of tables: %s\n' % len(tables))
#    for i, t in enumerate(tables):
#        #print(t)
#        #print(i)
#        #title = t.split('\n')[0][2:]
#        #t = '\n'.join(t.split('\n')[1:])
#        df = pd.read_csv(StringIO(t), sep=sep)
#        logfile.write(df, '\n\n')
#        stapp
#
#        # Parse table here, into oil object?
#        try:
#            parse_weathering_table(df, oil)
#        except:
#            logfile.write('Not a weathering table\n')
#
#    #print(tables, 'TAB')
#    logfile.write('Finished_parsing_tables\n')
#    jolla
#
#results = {}  # to store the parsed values
#z = zipfile.ZipFile(b+'.zip')
#logfile.write('Parsing zipfile ' + b+'.zip\n')
#for table in z.namelist():
#    #wd = {}
#    #print('Parsing table ' + table)
#    df = pd.read_csv(z.open(table))
#    if df.isnull().all().all():
#        continue  # empty table
#    logfile.write('# %s\n' % table)
#    logfile.write('%s\n\n' % df.to_string())
#    continue
#    if table != 'TOR_oil-page-95-table-1.csv':
#        continue
#
#    parse_weathering_table(df, oil)
#    parsing
#
#    # Search for row parameters
#    for v, r in template['rows'].items():
#        for name in r['names']:
#            for i, row in df.iterrows():
#                if name in str(row):
#                    varname = v
#                    while(varname) in results:
#                        varname = varname + 'a'
#                    results[varname] = []
#                    #print(varname + str(row.values))
#                    for ri in range(len(row.values)):
#                        try:
#                            row.values[ri] = np.float32(row.values[ri].replace(',', '.'))
#                        except:
#                            pass
#                    try:
#                        for aval in row.values:
#                            try:
#                                results[varname].append(np.float32(aval))
#                            except:
#                                pass
#                        results[varname] = np.array(results[varname])
#                        if 'factor' in r:
#                            results[varname] = results[varname]*r['factor']
#                        if 'subtract_from' in r:
#                            results[varname] = r['subtract_from'] - results[varname]
#                        if 'add_to' in r:
#                            results[varname] = r['add_to'] + results[varname]
#                        #print(results[varname])
#                    except:
#                        logfile.write('Could not parse row: %s\n' % str(row))
#                        stop
#                        pass
#    stop
#
#    #print('='*30)
#    #print(df)
#    #print('='*30)
#    for v, r in template['columns'].items():
#        headerline = df.columns
#        if 'Temp' in headerline[0] and ('Vol' in headerline[1] or 'Duva' in headerline[1]):
#            cut_temps = df.iloc[:,0].values + 273
#            cut_vols = df.iloc[:,1].values
#            try:
#                cut_vols = cut_vols / 100
#            except:
#                cut_vols = [np.float64(val.replace(',', '.'))/100 for val in cut_vols]
#
#            re = np.empty(len(cut_temps)*3)
#            for n in range(len(cut_temps)):
#                re[n*3] = cut_temps[n]
#                re[n*3+2] = cut_vols[n]
#            logfile.write('Cuts:\n' + np.array2string(re.astype(np.float32), separator='\t').replace('\n', '\t') + '\n')
##print('='*30)
##print(pd.read_csv(z.open(z.namelist()[0])))
##print('='*30)
#for n, r in results.items():
#    if len(r) > 0 and np.sum(np.isfinite(r))>0:
#        logfile.write(n, r, '\n')
#logfile.write('='*30 + '\n')
#
#for aas in ['', 'a', 'aa', 'aaa', 'aaaa', 'aaaaa']:
#    dstring = 'density' + aas
#    vstring = 'viscosity' + aas
#    wstring = 'weathered' + aas
#    rd = np.empty(12)
#    rv = np.empty(12)
#    for i in range(4):
#        if dstring in results and len(results[dstring])>0:
#            rd[i*3] = results[dstring][i]
#            rd[i*3+1] = 286
#            rd[i*3+2] = results[wstring][i]
#        if vstring in results and len(results[vstring])>0:
#            rv[i*3] = results[vstring][i]
#            rv[i*3+1] = 286
#            rv[i*3+2] = results[wstring][i]
#
#    if vstring in results:
#        logfile.write(vstring, np.array2string(rv.astype(np.float32), separator='\t').replace('\n', '\t')+ '\n')
#    if dstring in results:
#        logfile.write(dstring, np.array2string(rd.astype(np.float32), separator='\t').replace('\n', '\t') + '\n')
#
#
##try:
##    r = np.empty(12)
##    for i in range(4):
##        r[i*3] = results['density'][i]
##        r[i*3+1] = 286
##        r[i*3+2] = results['weathered'][i]
##    print(np.array2string(r.astype(np.float32), separator='\t').replace('\n', '\t'))
##except:
##    print('No density')
##
##try:
##    for i in range(4):
##        r[i*3] = results['viscosity'][i]
##    print(np.array2string(r.astype(np.float32), separator='\t').replace('\n', '\t'))
##except:
##    print('No viscosity')
