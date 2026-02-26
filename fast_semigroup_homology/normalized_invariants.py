"""
Given a finitely generated abelian group Z/c0 + Z/c1 + ... + Z/cn,
find an isomorphic abelian group Z/d0 + Z/d1 + ... + Z/dm
such that dm | ... | d1 | d0.

Note that we could have high numbers of some divisors here,
so we won't ever iterate through the entire list.
Instead, only operate on (divisor, count) pairs.
"""

from math import gcd
from itertools import pairwise

def _equivalent_invariant_factors(d1_count1, d2_count2):
    d1, count1 = d1_count1
    d2, count2 = d2_count2
    assert d1 > 1
    assert d2 > 1
    if d1 == d2:
        return [(d1, count1 + count2)]
    elif d2 % d1 == 0:
        return [d1_count1, d2_count2]
    elif d1 % d2 == 0:
        return [d2_count2, d1_count1]
    g = gcd(d1, d2)
    lcm = (d1//g)*d2
    if count1 < count2:
        if g == 1:
            return [(d2, count2-count1), (lcm, count1)]
        else:
            return [(g, count1), (d2, count2-count1), (lcm, count1)]
    elif count1 == count2:
        if g == 1:
            return [(lcm, count1)]
        else:
            return [(g, count1), (lcm, count1)]
    elif count1 > count2:
        if g == 1:
            return [(d1, count1-count2), (lcm, count2)]
        else:
            return [(g, count2), (d1, count1-count2), (lcm, count2)]

def invariant_factors(input_counter):
    free_rank = input_counter[0]
    data = [(d, count) for d, count in input_counter.items() if d > 1]
    data.sort()
    while not all(d1 < d2 and d2 % d1 == 0 for (d1, _), (d2, _) in pairwise(data)):
        i = 0
        while i < len(data) - 1:
            data[i:i+2] = _equivalent_invariant_factors(*data[i:i+2])
            i += 1
        data.sort()
    if free_rank:
        data.append((0, free_rank))
    return dict(reversed(data))
