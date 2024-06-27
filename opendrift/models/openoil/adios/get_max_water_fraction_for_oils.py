def get_water_fractions():
    import difflib
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from opendrift.models.openoil import OpenOil

    oot = OpenOil(location='Norway').oiltypes
    print(oot)

    o = pd.read_excel('oljedatabase_alle_03052024.xlsx')
    oljer = o['Oljetype'].unique().astype(list)
    oljer = oljer[0:-1]
    print(oljer)

    if False:
        for ot in oot:  # Loop to generate suggestions for mapping
            matches = difflib.get_close_matches(ot.title(), oljer, 1, cutoff=.5)
            if len(matches) > 0:
                matches = matches[0]
            else:
                matches = ''
            print(f'\'{ot}\': \'{matches}\',')

    mapping = {
        'AASGARD A 2003': 'Åsgard 2003',
        'AASTA HANSTEEN BLEND 2020': 'Aasta Hansteen Blend 2020',
        'ALTA 2016': 'Alta 2016',
        'ALVE 2010': 'Alve 2014 *(2017)',
        'ALVE 2014': 'Alve 2014 *(2017)',
        'ALVHEIM KNELER 2009': 'Kneler (Alvheim) 2009',
        'ATLA KONDENSAT 2013': 'Atla kondensat 2013',
        'AVALDSNES 2012': 'Avaldsnes (Johan Sverdrup) 2012',
        'BALDER 2002': 'Balder Blend 2010',
        'BALDER BLEND 2010': 'Balder Blend 2010',
        'BOYLA CRUDE 2016': 'Bøyla 2016',
        'BRAGE 2013': 'Brage 2013',
        'BRASSE 2018': 'Brasse 2018',
        'BREAM 2011': 'Bream 2011',
        'BREIDABLIKK 2023': 'Breidablikk',
        'BRYNHILD CRUDE 2015': 'Brynhild 2015',
        'CAURUS 2011': 'Caurus 2011',
        'DRAUGEN 2008': 'Draugen 2008',
        'DRIVIS 2017': 'Drivis 2017',
        'DUGONG 2022': 'Dugong 2022',
        'DUVA 2021': 'Duva 2021',
        'DVALIN 2020': 'Dvalin kondensat 2020',
        'EKOFISK 2002': 'Ekofisk blend 2002',
        'EKOFISK BLEND 2002': 'Ekofisk blend 2002',
        'EKOFISK BLEND 2015': 'Ekofisk J 2002',
        'EKOFISK J 2015': 'Ekofisk J 2015',
        'ELDFISK 2002': 'Eldfisk B 2015',
        'ELDFISK B 2015': 'Eldfisk B 2015',
        'ELDFISK BLEND 2015': 'Eldfisk B 2015',
        'ELDFISK KOMPLEKS 2015': 'Eldfisk Kompleks 2015',
        'ELLI 1999': 'Elli 1999',
        'ELLI SOUTH 1999': 'Elli South 1999',
        'EMBLA 2002': 'Embla 2002',
        'FENJA (PIL) 2015': 'Fenja (Pil) 2015',
        'FOGELBERG CONDENSATE 2021': 'Fogelberg Kondensat 2021',
        'FORSETI 2002': 'Forseti 2002',
        'FOSSEKALL 2013': 'Fossekall 2013',
        'FRAM 2013': 'Fram 2013',
        'FROSK 2020': 'Frosk 2020',
        'FROY 1996': 'Frøy 1996',
        'GARANTIANA 2013': 'Garantiana 2013',
        'GAUPE 2011': 'Gaupe 2011',
        'GINA KROG CRUDE 2018': 'Gina Krog 2018',
        'GJOA 2011': 'Gjøa 2011',
        'GLITNE 2002': 'Glitne 2002',
        'GOLIAT BLEND 50/50 2008': 'Goliat Blend 2008',
        'GOLIAT BLEND 70/30 2008': 'Goliat Blend 2008',
        'GOLIAT KOBBE 2008': 'Goliat Kobbe 2008 *(2010)',
        'GOLIAT REALGRUNNEN 2001': 'Goliat Realgrunnen 2003',
        'GOLIAT REALGRUNNEN 2008': 'Goliat Realgrunnen 2003',
        'GRANE 1997': 'Grane 1997',
        'GROSBEAK 2012': 'Grosbeak 2012',
        'GUDRUN 2012': 'Gudrun 2012',
        'GUDRUN 2019': 'Gudrun 2019',
        'GULLFAKS A BLEND 2010': 'Gullfaks A 2010',
        'GULLFAKS C BLEND 2010': 'Gullfaks C 2010',
        'GULLFAKS SOR 1996': 'Gullfaks Sør 1996',
        'GYDA 2002': 'Gyda 2002',
        'HAVIS 2013': 'Havis 2013',
        'HEIDRUN AaRE 2004': 'Heidrun Åre 2004',
        'HEIDRUN AARE 2023': 'Heidrun Åre 2023',
        'HEIDRUN EXPORT BLEND 2004': 'Heidrun Åre 2004',
        'HEIDRUN TILJE 2004': 'Heidrun Åre 2004',
        'HULDRA KONDENSAT 1998': 'Huldra 1998',
        'IRIS CONDENSATE 2020': 'Iris 2020',
        'IVAR AASEN 2012': 'Draupne (Ivar Aasen) 2012',
        'JORDBAER 2011': 'Jordbær 2011',
        'KRISTIN 2006': 'Kristin 2006',
        'KVITEBJORN 2009': 'Kvitebjørn 2009',
        'KVITEBJORN 2019': 'Kvitebjørn 2019',
        'LAVRANS 1997': 'Lavrans 1997',
        'LILLE PRINSEN 2022': 'Lille Prinsen 2022',
        'LILLEFRIGG KONDENSAT 1996': 'Lillefrigg 1996',
        'LINERLE 2005': 'Linerle 2005',
        'LUNO 2011': 'Luno (Edvard Grieg) 2011',
        'LUNO II 2014': 'Solveig (Luno II) 2014',
        'MARIA 2013': 'Maria 2013',
        'MARTIN LINGE CONDENSATE 2016': 'Martin Linge kondensat 2016',
        'MARTIN LINGE CRUDE 2016': 'Martin Linge olje 2016',
        'MARULK 2014': 'Marulk 2014',
        'MORVIN 2008': 'Morvin 2008',
        'NJORD 1997': 'Njord 1997* (2023)',
        'NJORD 2002': 'Njord 2003',
        'NJORD 2003': 'Njord 2003',
        'NORNE 2010': 'Norne 2018',
        'NORNE BLEND 2010': 'Norne Blend 2010',
        'NORNE CRUDE 2017': 'Norne 2018',
        'ODA 2019': 'Oda 2019',
        'OFELIA 2023': 'Ofelia 2023',
        'ORMEN LANGE KONDENSAT 2008': 'Ormen Lange 2008',
        'OSEBERG A 2013': 'Oseberg A 2013',
        'OSEBERG BLEND 2007': 'Oseberg Blend 2007',
        'OSEBERG C 1995': 'Oseberg C 1999',
        'OSEBERG C 2013': 'Oseberg A 2013',
        'OSEBERG OST 1998': 'Oseberg C 1999',
        'OSEBERG OST 2013': 'Oseberg Øst 2013',
        'OSEBERG SOR 2000': 'Oseberg Sør 2001',
        'OSEBERG SOR 2013': 'Oseberg Sør 2013',
        'OSELVAR 2012': 'Oselvar 2012',
        'RINGHORNE 2002': 'Ringhorne 2002',
        'SF NORD BRENT 2021': 'SF Nord Brent 2021',
        'SKARFJELL 2014': 'Skarfjell 2014',
        'SKARV 2004': 'Skarv olje 2004',
        'SKARV KONDENSAT 2014': 'Skarv kondensat 2014',
        'SKOGUL 2020': 'Skogul 2020',
        'SKRUGARD 2012': 'Skrugard 2012',
        'SLEIPNER KONDENSAT 2002': 'leipner 2002',
        'SLEIPNER VEST 1998': 'Sleipner Vest 1998',
        'SMORBUKK 2003': 'Smørbukk Sør 2003',
        'SMORBUKK KONDENSAT 2003': 'Smørbukk kondensat 2003',
        'SMORBUKK SOR 2003': 'Smørbukk Sør 2003',
        'SNOHVIT KONDENSAT 2001': 'Snøhvit 2001',
        'SNORRE B 2004': 'Snorre B 2004',
        'SNORRE TLP 2004': 'Snorre TLP 2004',
        'STAER 2010': 'Stær 2010',
        'STATFJORD A 2001': 'Statfjord A 2001',
        'STATFJORD B 2001': 'Statfjord A 2001',
        'STATFJORD C 2001': 'Statfjord C 2001 *(2021)',
        'SVALE 2010': 'Svale 2010',
        'SVALIN 2014': 'Svalin 2014',
        'SYGNA BRENT 2021': 'Sygna Brent 2021',
        'TAMBAR 2002': 'Tambar 2002',
        'TAU 1999': 'Tau 1999',
        'TOR 2002': 'Tambar 2002',
        'TOR II 2022': 'TOR ll 2022',
        'TORDIS 2002': 'Iris 2020',
        'TRESTAKK 2008': 'Trestakk 2008',
        'TROLL, STATOIL': 'Troll 1999',
        'TRYM KONDENSAT 2011': 'Trym 2011',
        'TYRIHANS NORD 2004': 'Tyrihans Nord 2004',
        'TYRIHANS SOR 2004': 'Tyrihans Sør 2004',
        'ULA 1999': 'Ula 1999',
        'UTGARD CONDENSATE 2021': 'Utgard kondensat 2021',
        'VALE 2001': 'Vale 2001',
        'VALE 2014': 'Vale 2014',
        'VALHALL 2002': 'Valhall 2002',
        'VALHALL 2021': 'Valhall 2021',
        'VARG 2000': 'Varg 2000',
        'VEGA CONDENSATE 2015': 'Vega 2016',
        'VESLEFRIKK 2012': 'Veslefrikk 2012',
        'VILJE 2009': 'Vilje 2009',
        'VISUND 2009': 'Visund 2009',
        'VISUND CRUDE OIL 2020': 'Visund Crude Oil 2020',
        'VISUND SOR CONDENSATE 2020': 'Visund Sør Kondensat 2020',
        'VOLUND 2010': 'Volund 2010',
        'VOLVE 2006': 'Volve 2006',
        'WISTING 2015': 'Wisting 2015',
        'WISTING CENTRAL 2017': 'Wisting Central 2017',
        'YME 2023': 'Yme 2023'
        }


    #o = o.loc[o['Oljetype'] == 'Yme 2023']
    print(o)
    print(o.columns)
    print(o.groupby(['Oljetype', 'Temperatur (°C)'])['Vanninnhold'].max().to_string())

    oilmax = {  # Some known limits not in NOFO excel file
            'MARINE GAS OIL 500 ppm S 2017':
                {'temperatures': [15], 'max_water_fraction': [.1]},
            }

    for ot in mapping:
        d = o.loc[o['Oljetype'] == mapping[ot]].groupby(['Temperatur (°C)'])['Vanninnhold'].max()
        oilmax[ot] = {'temperatures': list(d.keys().values),
                      'max_water_fraction': list(d.values)}
    print(oilmax)
    with open('max_water_fraction.json', 'w') as json_file:
       json.dump(oilmax , json_file)

if __name__ == '__main__':
    get_water_fractions()
