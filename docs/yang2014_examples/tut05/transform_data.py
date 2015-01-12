from __future__ import print_function, division, absolute_import

from collections import defaultdict
import string


def canonical(seq_in):
    d = {}
    seq_out = []
    for x in seq_in:
        if x not in d:
            d[x] = len(d)
        seq_out.append(d[x])
    return seq_out


def get_info():
    aas = ''.join(x for x in string.ascii_uppercase if x not in 'BJOUXZ')
    d = {a : i for i, a in enumerate(aas)}

    sequences = []
    with open('mtCDNApri.aa') as fin:
        lines = fin.readlines()
        header = lines[0]
        for line in lines[1:]:
            name, sequence = line.strip().split()
            sequences.append([d[x] for x in sequence])

    canonical_columns = defaultdict(int)
    for column in zip(*sequences):
        ccol = canonical(column)
        canonical_columns[tuple(ccol)] += 1

    return sequences, canonical_columns


def main():
    sequences, canonical_columns = get_info()

    lines = []
    for column in zip(*sequences):
        lines.append(str(list(column)))

    print('all columns:')
    print(',\n'.join(lines))
    print()

    canonical_lines = []
    counts = []
    pairs = sorted(canonical_columns.items())
    for ccol, count in pairs:
        canonical_lines.append(str(list(ccol)))
        counts.append(count)

    print('canonical columns and corresponding counts:')
    print(',\n'.join(canonical_lines))
    print(counts)
    print(range(len(counts)))
    print()

    print(len(lines), 'total columns')
    print(len(canonical_columns), 'columns unique up to residue label')


if __name__ == '__main__':
    main()
