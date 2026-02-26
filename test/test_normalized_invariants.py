from collections import Counter
from fast_semigroup_homology.normalized_invariants import invariant_factors
import random
from itertools import pairwise
from math import gcd

import pytest

def invariant_factors_simple_list(divisors):
    free_rank = divisors.count(0)
    divisors = list(filter(None, divisors))
    divisors.sort()
    while not all(d2 % d1 == 0 for d1, d2 in pairwise(divisors)):
        for i in range(len(divisors) - 1):
            d1, d2 = divisors[i:i+2]
            g = gcd(d1, d2)
            lcm = (d1//g)*d2
            divisors[i:i+2] = [g, lcm]
        divisors.sort()
    return [d for d in divisors if d > 1] + [0] * free_rank

@pytest.mark.parametrize("input,expected", [
    ({}, {}),
    ({0: 1}, {0: 1}),
    ({2: 7, 3: 10}, {6: 7, 3: 3}),
    ({4: 7, 10: 10}, {20: 7, 10: 3, 2: 7}),
    ({2: 10, 3: 10}, {6: 10}),
    ({4: 10, 10: 10}, {20: 10, 2: 10}),
    ({2: 10, 3: 7}, {6: 7, 2: 3}),
    ({4: 10, 10: 7}, {20: 7, 4: 3, 2: 7}),
    ({2: 3, 4: 1, 3: 10, 0: 50}, {0: 50, 12: 1, 6: 3, 3: 6}),
])
def test_invariant_factors(input, expected):
    input = Counter(input)
    result = invariant_factors(input)
    assert result == expected
    assert result == dict(Counter(invariant_factors_simple_list(list(input.elements()))))
    # order should also be preserved
    assert list(result.items()) == list(expected.items())

    N = 1_000_000 # Large number should do all the same things
    big_input = Counter({d: N*count for d, count in input.items()})
    big_expected = {d: N*count for d, count in expected.items()}
    big_result = invariant_factors(big_input)
    assert big_result == big_expected
    assert list(big_result.items()) == list(big_expected.items())

def test_invariant_factors_random_small_numbers():
    for _ in range(10):
        data0 = [random.getrandbits(4) for _ in range(30)]
        result = invariant_factors(Counter(data0))
        expected = dict(Counter(reversed(invariant_factors_simple_list(data0))))
        assert result == expected, data0
        assert list(result.items()) == list(expected.items()), data0
