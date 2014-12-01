from __future__ import print_function

from math import exp

def p0(mu, lam, t):
    a = exp(-(lam - mu)*t)
    return (mu * (1-a)) / (lam - mu*a)

def p1(mu, lam, t):
    a = exp(-(lam - mu)*t)
    u = (lam - mu) * (lam - mu) * a
    v = lam - mu * a
    return u / (v*v)

def pn(mu, lam, t, n):
    if n == 0:
        return p0(mu, lam, t)
    if n == 1:
        return p1(mu, lam, t)
    return (lam / mu)**(n-1) * p1(mu, lam, t) * p0(mu, lam, t)**(n-1)

def main():
    n = 1000
    mu = 0.3
    lam = 0.5
    t = 5
    print('mu:', mu)
    print('lam:', lam)
    expectation = 0
    p_failed = pn(mu, lam, t, 0)
    p_counted = 0
    for i in range(1, n):
        p = pn(mu, lam, t, i)
        p_counted += p
        expectation += i * p
    expectation /= p_counted
    print('probability of no survival:', p_failed)
    print('probability of reasonable population:', p_counted)
    print('probability of huge population:', 1 - (p_failed + p_counted))
    print('posterior expectation conditional on survival:', expectation)

main()
