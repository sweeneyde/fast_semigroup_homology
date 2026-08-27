import random
import itertools

from mutable_lattice import Vector, Lattice, transpose, relations_among
from math import gcd
from fast_semigroup_homology.cup_products.homology_with_generators import (
    _smith_with_coefficients_unnormalized,
    _normalize_pair,
    _normalize_smith_with_coefficients,
    _smith_with_coefficients,
    _smith_with_right_coefficients,
    _inverse,
    Cokernel,
    homology,
)

def test__inverse():
    assert _inverse([]) == []
    assert _inverse([Vector([1])]) == [Vector([1])]
    assert _inverse([Vector([-1])]) == [Vector([-1])]
    assert _inverse([Vector([1,0]),Vector([0,1])]) == [Vector([1,0]),Vector([0,1])]
    assert _inverse([Vector([1,0]),Vector([0,-1])]) == [Vector([1,0]),Vector([0,-1])]
    assert _inverse([Vector([1,5]),Vector([0,1])]) == [Vector([1,-5]),Vector([0,1])]
    assert _inverse([Vector([1,5]),Vector([0,1])]) == [Vector([1,-5]),Vector([0,1])]
    assert _inverse([Vector([9,4]),Vector([11,5])]) == [Vector([5,-4]),Vector([-11,9])]
    assert _inverse([Vector([11,5]),Vector([9,4])]) == [Vector([-4,5]),Vector([9,-11])]

def test_smith_with_coefficients_unnormalized():
    assert _smith_with_coefficients_unnormalized(0, []) == ([], [])
    assert _smith_with_coefficients_unnormalized(1, []) == ([], [])
    assert _smith_with_coefficients_unnormalized(2, []) == ([], [])
    assert _smith_with_coefficients_unnormalized(1, [Vector([1])]) == ([1], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(1, [Vector([-1])]) == ([1], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(1, [Vector([5])]) == ([5], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(1, [Vector([-5])]) == ([5], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(1, [Vector([0])]) == ([0], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(1, [Vector([5]), Vector([5])]) == ([5, 0], [Vector([0,1]),Vector([1,-1])])
    assert _smith_with_coefficients_unnormalized(2, [Vector([0,0])]) == ([0], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(2, [Vector([1,0])]) == ([1], [Vector([1])])
    assert _smith_with_coefficients_unnormalized(2, [Vector([1,0]), Vector([2, 2])]) == ([1, 2], [Vector([1,0]),Vector([0,1])])
    assert _smith_with_coefficients_unnormalized(2, [Vector([1,0]), Vector([2, 2]), Vector([0,2])]) == ([1, 2, 0], [Vector([1,0,0]),Vector([0,1,0]),Vector([2,-1,1])])
    assert _smith_with_coefficients_unnormalized(2, [Vector([2,0]), Vector([3, 3])]) == ([2, 3], [Vector([1,0]),Vector([0,1])])
    assert _smith_with_coefficients_unnormalized(2, [Vector([-5,0]), Vector([5,4])]) == ([1, 20], [Vector([-1,2]),Vector([-3,5])])

def test_normalize_pair():
    assert _normalize_pair(2, Vector([10,0]), 3, Vector([0,10])) == (1, Vector([-10,10]), 6, Vector([-30,20]))

def test_normalize_pair_random():
    vectors = [Vector([x, x+2]) for x in range(-3, 4)]
    pairs = [(d1, d2) for d1 in range(1, 7) for d2 in range(1, 7) if d2 % d1 != 0]
    for v1 in vectors:
        for v2 in vectors:
            for d1, d2 in pairs:
                d3, v3, d4, v4 = _normalize_pair(d1, v1, d2, v2)
                assert Lattice(2, [v1, v2]) == Lattice(2, [v3, v4])
                assert d3 == gcd(d1, d2)
                assert d4 == d1*d2 // d3

def test__smith_with_coefficients():
    assert _smith_with_coefficients(0, []) == ([], [])
    assert _smith_with_coefficients(1, []) == ([], [])
    assert _smith_with_coefficients(2, []) == ([], [])
    assert _smith_with_coefficients(1, [Vector([1])]) == ([1], [Vector([1])])
    assert _smith_with_coefficients(1, [Vector([-1])]) == ([1], [Vector([1])])
    assert _smith_with_coefficients(1, [Vector([5])]) == ([5], [Vector([1])])
    assert _smith_with_coefficients(1, [Vector([-5])]) == ([5], [Vector([1])])
    assert _smith_with_coefficients(1, [Vector([0])]) == ([0], [Vector([1])])
    assert _smith_with_coefficients(1, [Vector([5]), Vector([5])]) == ([5, 0], [Vector([0,1]),Vector([1,-1])])
    assert _smith_with_coefficients(2, [Vector([0,0])]) == ([0], [Vector([1])])
    assert _smith_with_coefficients(2, [Vector([1,0])]) == ([1], [Vector([1])])
    assert _smith_with_coefficients(2, [Vector([1,0]), Vector([2,2])]) == ([1, 2], [Vector([1,0]),Vector([0,1])])
    assert _smith_with_coefficients(2, [Vector([1,0]), Vector([2,2]), Vector([0,2])]) == ([1, 2, 0], [Vector([1,0,0]),Vector([0,1,0]),Vector([2,-1,1])])
    assert _smith_with_coefficients(2, [Vector([2,0]), Vector([3,3])]) == ([1, 6], [Vector([-1,1]),Vector([-3,2])])
    assert _smith_with_coefficients(2, [Vector([-5,0]), Vector([5,4])]) == ([1, 20], [Vector([-1,2]),Vector([-3,5])])

def test_smith_with_coefficients_random():
    for N, R, _ in itertools.product(range(5), repeat=3):
        data = [Vector([random.randint(-10, 10) for _ in range(N)]) for _ in range(R)]
        invariants, coefficient_vectors = _smith_with_coefficients(N, data)
        nonzero_invariants = [d for d in invariants if d]
        assert invariants == Lattice(R, transpose(N, data)).invariants()
        assert nonzero_invariants == Lattice(N, data).nonzero_invariants()
        # transformation was invertible
        assert Lattice(R, coefficient_vectors).invariants() == [1] * R
        combos = []
        for d, w in zip(invariants, coefficient_vectors):
            v = Vector.zero(N)
            for c, vi in zip(w, data):
                v += c * vi
            if d == 0:
                assert not any(v)
            else:
                assert Lattice(N, [v], maxrank=1).nonzero_invariants() == [d]
                combos.append(v)
        assert Lattice(N, combos).nonzero_invariants() == nonzero_invariants

def test_cokernel_examples():
    c = Cokernel(0, [])
    assert c.get_generators() == []
    assert c.get_invariants() == []
    c = Cokernel(0, [Vector([])])
    assert c.get_generators() == []
    assert c.get_invariants() == []
    c = Cokernel(0, [Vector([]), Vector([])])
    assert c.get_generators() == []
    assert c.get_invariants() == []

    c = Cokernel(1, [])
    assert c.get_generators() == [Vector([1])]
    assert c.get_invariants() == [0]
    c = Cokernel(1, [Vector([0])])
    assert c.get_generators() == [Vector([1])]
    assert c.get_invariants() == [0]
    c = Cokernel(1, [Vector([1])])
    assert c.get_generators() == []
    assert c.get_invariants() == []
    c = Cokernel(1, [Vector([-1])])
    assert c.get_generators() == []
    assert c.get_invariants() == []
    c = Cokernel(1, [Vector([5])])
    assert c.get_generators() == [Vector([1])]
    assert c.get_invariants() == [5]
    c = Cokernel(1, [Vector([-5])])
    assert c.get_generators() == [Vector([1])]
    assert c.get_invariants() == [5]

    c = Cokernel(2, [])
    assert c.get_generators() == [Vector([1, 0]), Vector([0, 1])]
    assert c.get_invariants() == [0, 0]
    c = Cokernel(2, [Vector([0, 0])])
    assert c.get_generators() == [Vector([1, 0]), Vector([0, 1])]
    assert c.get_invariants() == [0, 0]
    c = Cokernel(2, [Vector([1, 0])])
    assert c.get_generators() == [Vector([0, 1])]
    assert c.get_invariants() == [0]
    c = Cokernel(2, [Vector([2, -1])])
    assert c.get_generators() == [Vector([1, 0])]
    assert c.get_invariants() == [0]
    c = Cokernel(2, [Vector([1, 0]),Vector([2, 2])])
    assert c.get_generators() == [Vector([0, 1])]
    assert c.get_invariants() == [2]
    c = Cokernel(2, [Vector([1, 0]),Vector([2, 2]),Vector([0, 2])])
    assert c.get_generators() == [Vector([0, 1])]
    assert c.get_invariants() == [2]
    c = Cokernel(2, [Vector([2, 0]),Vector([3, 3])])
    assert c.get_generators() == [Vector([0, 1])]
    assert c.get_invariants() == [6]
    c = Cokernel(2, [Vector([-5, 0]),Vector([5, 4])])
    assert c.get_generators() == [Vector([-1, -1])]
    assert c.get_invariants() == [20]

    c = Cokernel(3, [Vector([1,0,0]), Vector([0,0,1])])
    assert c.get_invariants() == [0]
    assert c.get_generators() == [Vector([0, 1, 0])]
    assert c.projection(Vector([123, 456, 789])) == [456]

    c = Cokernel(3, [Vector([1,0,0]), Vector([0, 100, 0]), Vector([0,0,1])])
    assert c.get_invariants() == [100]
    assert c.get_generators() == [Vector([0, 1, 0])]
    assert c.projection(Vector([123, 456, 789])) == [56]

    c = Cokernel(3, [Vector([0,2,2]), Vector([0,0,4])])
    assert c.get_invariants() == [2, 4, 0]
    assert c.get_generators() == [Vector([0, 1, 1]), Vector([0, 0, 1]), Vector([1, 0, 0])]
    assert c.projection(Vector([2, 3, 5])) == [1, 2, 2]
    assert c.projection(Vector([7, 7, 7])) == [1, 0, 7]

def check_exponents(L, invariants, generators):
    BIG = 2520  # lcm(*range(1,11))

    # check that the invariant is the smallest coefficient for the generator
    # that makes it land in L
    for d, v in zip(invariants, generators):
        for k in range(1, d):
            assert k*v not in L
        assert d*v in L
        if d == 0:
            assert BIG*v not in L

    # Check for independence:
    # adding in the other vectors shouldn't change the exponent
    L1 = L.copy()
    for d, v in zip(invariants, generators):
        for k in range(1, d):
            assert k*v not in L1
        assert d*v in L1
        if d == 0:
             assert BIG*v not in L1
        L1.add_vector(v)

    # Check in the reverse direction also
    L2 = L.copy()
    for d, v in zip(invariants[::-1], generators[::-1]):
        for k in range(1, d):
            assert k*v not in L2
        assert d*v in L2
        if d == 0:
             assert BIG*v not in L2
        L2.add_vector(v)

def test_cokernel_random():
    for N, R, _ in itertools.product(range(5), repeat=3):
        data = [Vector([random.randint(-10, 10) for _ in range(N)]) for _ in range(R)]
        c = Cokernel(N, data)
        generators = c.get_generators()
        invariants = c.get_invariants()
        L = Lattice(N, data)
        assert c.get_invariants() == [d for d in L.invariants() if d != 1]
        assert (L + Lattice(N, generators)).is_full()
        check_exponents(L, invariants, generators)
        # projection after inclusion (retraction)
        for i, gen in enumerate(generators):
            assert c.projection(gen) == [0] * i + [1] + [0] * (len(generators)-i-1), data
        # inclusion after projection agrees up to L elements
        zero = Vector.zero(N)
        for v in data:
            v1 = sum((x*gen for x, gen in zip(c.projection(v), generators)), start=zero)
            assert v1 - v in L

def test_homology_examples():
    h = homology(
        [Vector([0,0]), Vector([0,0]), Vector([0,0])],
        [Vector([7,7,7,7]),Vector([0,0,5,5])]
    )
    assert h.get_generators() == []
    assert h.get_invariants() == []
    h = homology(
        [Vector([0,0])],
        [Vector([7,7,7,7]),Vector([7,7,7,7])]
    )
    assert h.get_generators() == [Vector([1,-1])]
    assert h.get_invariants() == [0]
    h = homology(
        [Vector([0,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3])]
    )
    assert h.get_generators() == [Vector([3,-2])]
    assert h.get_invariants() == [0]
    h = homology(
        [Vector([0,0,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_generators() == [Vector([3,-2, 0]), Vector([0, 0, 1])]
    assert h.get_invariants() == [0, 0]
    h = homology(
        [Vector([3,-2,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_generators() == [Vector([0, 0, 1])]
    assert h.get_invariants() == [0]
    h = homology(
        [Vector([0,0,1])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_generators() == [Vector([3,-2,0])]
    assert h.get_invariants() == [0]
    h = homology(
        [Vector([0,0,1]), Vector([3,-2,1])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_generators() == []
    assert h.get_invariants() == []
    h = homology(
        [Vector([3,-2,-100])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_generators() == [Vector([0,0,1])]
    assert h.get_invariants() == [0]
    h = homology(
        [Vector([-30,20,1])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_generators() == [Vector([3,-2,0])]
    assert h.get_invariants() == [0]

    h = homology(
        [Vector([0,0])],
        [Vector([7,7,7,7]),Vector([7,7,7,7])]
    )
    assert h.get_invariants() == [0]
    assert h.get_generators() == [Vector([1, -1])]
    assert h.projection(Vector([2,-2])) == [2]

    h = homology(
        [Vector([3,-2,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert h.get_invariants() == [0]
    assert h.get_generators() == [Vector([0, 0, 1])]
    assert h.projection(Vector([0, 0, 1])) == [1]
    assert h.projection(Vector([3, -2, 1])) == [1]
    assert h.projection(Vector([300, -200, 17])) == [17]


def test_homology_random():
    for K, N, R, _ in itertools.product(range(5), repeat=4):
        outgoing = [Vector([random.randint(-10,10) for _ in range(K)])
                    for _ in range(N)]
        kernel = relations_among(outgoing)
        incoming = [kernel.linear_combination(
                        Vector([random.randint(-10,10)
                                for _ in range(kernel.rank)]))
                    for _ in range(R)]
        h = homology(incoming, outgoing)
        generators = h.get_generators()
        invariants = h.get_invariants()
        L = Lattice(N, incoming)
        assert h.get_invariants() == [d for d in L.nonzero_invariants() if d != 1] + [0] * (kernel.rank - L.rank)
        assert L + Lattice(N, generators) == kernel
        check_exponents(L, invariants, generators)
        # projection after inclusion (retraction)
        for i, gen in enumerate(generators):
            assert h.projection(gen) == [0] * i + [1] + [0] * (len(generators)-i-1)
        # inclusion after projection agrees up to L elements
        zero = Vector.zero(N)
        for v in incoming:
            v1 = sum((x*gen for x, gen in zip(h.projection(v), generators)), start=zero)
            assert v1 - v in L
