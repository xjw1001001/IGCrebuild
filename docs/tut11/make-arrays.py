from __future__ import print_function, division

import itertools

def get_diffs(a, b):
    return [(x, y) for x, y in zip(a, b) if x != y]

def main():

    # List all of the self-consistent multivariate states.
    pri_to_tol = [0, 0, 1, 1, 2, 2]
    legit_states = []
    for state in itertools.product([0, 1, 2, 3, 4, 5], [0, 1], [0, 1], [0, 1]):
        pri, tols = state[0], state[1:]
        if not tols[pri_to_tol[pri]]:
            continue
        s = list(state)
        legit_states.append(s)

    # print all legit states
    for s in legit_states:
        print(s)
    print()

    # List all of the possible transitions.
    forward_trans = [
            (0, 1), (2, 3), (4, 5),
            (0, 2), (2, 4),
            (1, 3), (3, 5),
            ]
    backward_trans = [(b, a) for a, b in forward_trans]
    trans = forward_trans + backward_trans
    rows = []
    cols = []
    rates = []
    for sa in legit_states:
        for sb in legit_states:
            diffs = get_diffs(sa, sb)
            if len(diffs) != 1:
                continue
            d = diffs[0]
            if d in trans:
                rows.append(sa)
                cols.append(sb)
                rates.append(1)

    print(',\n'.join(str(s) for s in rows))
    print()

    print(',\n'.join(str(s) for s in cols))
    print()

    print(rates)

main()
