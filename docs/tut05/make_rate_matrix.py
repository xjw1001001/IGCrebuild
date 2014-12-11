from __future__ import print_function, division

def gen_triples(n):
    for i in range(n):
        for j in range(n):
            if i != j:
                yield i, j, 1

def main(n):

    n = 20

    print('row_states')
    lines = []
    for i in range(n):
        line = ','.join(('[%d]' % i) for j in range(n) if i != j)
        lines.append(line)
    print(',\n'.join(lines))

    print('column_states')
    lines = []
    for i in range(n):
        line = ','.join(('[%d]' % j) for j in range(n) if i != j)
        lines.append(line)
    print(',\n'.join(lines))

    print('transition_rates')
    lines = []
    for i in range(n):
        line = ','.join('1' for j in range(n) if i != j)
        lines.append(line)
    print(',\n'.join(lines))


if __name__ == '__main__':
    main()
