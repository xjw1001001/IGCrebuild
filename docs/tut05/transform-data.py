from __future__ import print_function, division

import string

aas = ''.join(x for x in string.ascii_uppercase if x not in 'BJOUXZ')
d = {a : i for i, a in enumerate(aas)}

sequences = []
with open('mtCDNApri.aa') as fin:
    lines = fin.readlines()
    header = lines[0]
    for line in lines[1:]:
        name, sequence = line.strip().split()
        sequences.append([d[x] for x in sequence])

lines = []
for column in zip(*sequences):
    lines.append(str(list(column)))
print(',\n'.join(lines))
