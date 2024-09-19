#!/usr/bin/env python

import numpy as np

f = open('/home/knutfd/software/OilLibrary/oil_library/OilLibNorway', 'r')
o = open('oil_types.rst', 'w')

o.write('Oil types\n#########\n')
o.write('In the table below, some key parameters and reference is given for the Norwegian oil types of the NOAA OilLibrary.\n\n')

o.write('.. list-table::\n   :widths: 20 10 10 10 10 30\n   :header-rows: 1\n\n')
o.write('   * - Name\n     - Density [kg/m3]\n     - Viscosity [kg*m-1*s-1]\n     - Asphaltenes [fraction]\n     - Wax cont [fraction]\n     - Reference\n\n')

for line in f.readlines()[2::]:
    el = line.split('\t')
    if el[23] == '':
        den = np.nan
    else:
        den = float(el[23])
    name = el[0]
    if el[53] == '':
        vis = np.nan
    else:
        vis = float(el[53])
    asp = el[11]
    wax = el[12]
    ref = el[5]
    
    o.write('   * - %s\n     - %s\n     -  %s\n     -  %s\n     -  %s\n     -  %s\n\n' % (name, den, vis, asp, wax, ref))
