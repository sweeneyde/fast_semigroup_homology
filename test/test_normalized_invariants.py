from collections import Counter
from fast_monoid_homology.normalized_invariants import invariant_factors
import random
from itertools import pairwise
from math import gcd

def test_invariant_factors():
    assert invariant_factors(Counter({})) == {}
    assert invariant_factors(Counter({0: 1})) == {0: 1}
    assert invariant_factors(Counter({2: 1000})) == {2: 1000}
    assert invariant_factors(Counter({2: 3_000, 3: 10_000})) == {3: 7_000, 6: 3_000}
    assert invariant_factors(Counter({2: 3_000, 4: 1_000, 3: 10_000})) == {3: 6_000, 6: 3_000, 12: 1_000}

def test_invariant_factors_random_small_numbers():
    data0 = [random.getrandbits(8) for _ in range(300)]
    result = invariant_factors(Counter(data0))
    data = data0[:]
    data.sort()
    data.sort(key=lambda x: x==0)
    while not all(d2 == 0 or d2 % d1 == 0 for d1, d2 in pairwise(data)):
        for i in range(len(data) - 1):
            d1, d2 = data[i:i+2]
            if d2 != 0 and d2 % d1 != 0:
                g = gcd(d1, d2)
                lcm = (d1//g)*d2
                data[i:i+2] = g, lcm
        data = [x for x in data if x != 1]
        data.sort()
        data.sort(key=lambda x: x==0)
    assert dict(Counter(data)) == result, data0
