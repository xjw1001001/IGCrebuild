def get_MG94BasicRate(ca, cb, pi, kappa, omega, codon_table):
    dif = [ii for ii in range(3) if ca[ii] != cb[ii]]
    ndiff = len(dif)
    if ndiff > 1:
        return 0
    elif ndiff == 0:
        print 'Please check your codon tables and make sure no redundancy'
        print ca, cb
        return 0
    else:
        na = ca[dif[0]]
        nb = cb[dif[0]]
        QbasicRate = pi['ACGT'.index(nb)]

        if isTransition(na, nb):
            QbasicRate *= kappa

        if isNonsynonymous(ca, cb, codon_table):
            QbasicRate *= omega

        return QbasicRate

def isTransition(na, nb):
    return (set([na, nb]) == set(['A', 'G']) or set([na, nb]) == set(['C', 'T']))

def isNonsynonymous(ca, cb, codon_table):
    return (codon_table[ca] != codon_table[cb])

vec_get_MG94BasicRate = np.vectorize(get_MG94BasicRate, doc='Vectorized `get_MG94BasicRate`', excluded = ['pi', 'kappa', 'omega', 'codon_table'])
