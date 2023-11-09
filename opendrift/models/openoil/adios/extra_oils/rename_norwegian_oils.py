import json
import glob
import adios_db.scripting as ads


files = glob.glob('AD03*.json')
print(files)

for f in files:
    o = json.load(open(f))
    o['oil_id'] = o['oil_id'].replace('AD03', 'NO00')
    f = f.replace('AD03', 'NO00')
    o['metadata']['source_id'] = o['metadata']['source_id'].replace('AD03', 'NO00')
    name = o['metadata']['name']
    year = o['metadata']['reference']['year']
    if name[-1].isdigit() and name[-4].isdigit():
        newname = name[0:-5].strip()
    elif name[-1] == ')' and name[-6] == '(':
        newname = name[0:-6].strip()
    else:
        newname = 'unchanged'
    o['metadata']['name'] = newname

    print(f'*{name}* - *{newname}* - {year}')
    print('\n\n\n')

    no = ads.Oil.from_py_json(o)
    no.to_file(f)

