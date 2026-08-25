import random
import itertools

from mutable_lattice import Vector, Lattice, transpose, relations_among
from math import gcd
from fast_semigroup_homology.cup_products.homology_with_generators import (
    _smith_with_coefficients_unnormalized,
    _normalize_pair,
    _normalize_smith_with_coefficients,
    _smith_with_coefficients,
    _cokernel_with_generators_unfiltered,
    cokernel_with_generators,
    _inverse,
    homology_with_generators,
    cokernel_with_generators_and_projection,
    homology_with_generators_and_projection,
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

def test_cokernel_with_generators():
    assert cokernel_with_generators(0, []) == ([], [])
    assert cokernel_with_generators(1, []) == ([0], [Vector([1])])
    assert cokernel_with_generators(2, []) == ([0, 0], [Vector([1,0]), Vector([0,1])])
    assert cokernel_with_generators(1, [Vector([1])]) == ([], [])
    assert cokernel_with_generators(1, [Vector([-1])]) == ([], [])
    assert cokernel_with_generators(1, [Vector([5])]) == ([5], [Vector([1])])
    assert cokernel_with_generators(1, [Vector([-5])]) == ([5], [Vector([1])])
    assert cokernel_with_generators(1, [Vector([0])]) == ([0], [Vector([1])])
    assert cokernel_with_generators(1, [Vector([5]), Vector([5])]) == ([5], [Vector([1])])
    assert cokernel_with_generators(2, [Vector([0,0])]) == ([0, 0], [Vector([1,0]),Vector([0,1])])
    assert cokernel_with_generators(2, [Vector([1,0])]) == ([0], [Vector([0,1])])
    assert cokernel_with_generators(2, [Vector([2,-1])]) == ([0], [Vector([1,0])])
    assert cokernel_with_generators(2, [Vector([1,0]),Vector([2,2])]) == ([2], [Vector([0,1])])
    assert cokernel_with_generators(2, [Vector([1,0]),Vector([2,2]),Vector([0,2])]) == ([2], [Vector([0,1])])
    assert cokernel_with_generators(2, [Vector([2,0]),Vector([3,3])]) == ([6], [Vector([0,1])])
    assert cokernel_with_generators(2, [Vector([-5,0]),Vector([5,4])]) == ([20], [Vector([-1,-1])])

def test_cokernel_with_generators_random():
    for N, R, _ in itertools.product(range(5), repeat=3):
        data = [Vector([random.randint(-10, 10) for _ in range(N)]) for _ in range(R)]
        invariants, generators = _cokernel_with_generators_unfiltered(N, data)
        L = Lattice(N, data)
        L0 = L.copy()
        assert invariants == L.nonzero_invariants() + [0] * (N - L.rank)
        for d, v in zip(invariants, generators):
            for k in range(1, d):
                assert k*v not in L
                assert k*v not in L0
            assert d*v in L
            assert d*v in L0
            L.add_vector(v)
        assert L.is_full()

def test_homology_with_generators():
    assert homology_with_generators(
        [Vector([0,0]), Vector([0,0]), Vector([0,0])],
        [Vector([7,7,7,7]),Vector([0,0,5,5])]
    ) == ([], [])
    assert homology_with_generators(
        [Vector([0,0])],
        [Vector([7,7,7,7]),Vector([7,7,7,7])]
    ) == ([0], [Vector([1,-1])])
    assert homology_with_generators(
        [Vector([0,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3])]
    ) == ([0], [Vector([3,-2])])
    assert homology_with_generators(
        [Vector([0,0,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    ) == ([0, 0], [Vector([3,-2, 0]), Vector([0, 0, 1])])
    assert homology_with_generators(
        [Vector([3,-2,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    ) == ([0], [Vector([0, 0, 1])])
    assert homology_with_generators(
        [Vector([0,0,1])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    ) == ([0], [Vector([3,-2,0])])
    assert homology_with_generators(
        [Vector([0,0,1]), Vector([3,-2,1])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    ) == ([], [])
    assert homology_with_generators(
        [Vector([3,-2,-100])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    ) == ([0], [Vector([0,0,1])])
    assert homology_with_generators(
        [Vector([-30,20,1])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    ) == ([0], [Vector([3,-2,0])])

def test_homology_with_generators_random():
    for K, N, R, _ in itertools.product(range(5), repeat=4):
        outgoing = [Vector([random.randint(-10,10) for _ in range(K)])
                    for _ in range(N)]
        kernel = relations_among(outgoing)
        incoming = [kernel.linear_combination(
                        Vector([random.randint(-10,10)
                                for _ in range(kernel.rank)]))
                    for _ in range(R)]
        invariants, generators = homology_with_generators(incoming, outgoing)
        L = Lattice(N, incoming)
        L0 = L.copy()
        assert invariants == [d for d in L.nonzero_invariants() if d != 1] + [0] * (kernel.rank - L.rank)
        for d, v in zip(invariants, generators):
            for k in range(1, d):
                assert k*v not in L
                assert k*v not in L0
            assert d*v in L
            assert d*v in L0
            L.add_vector(v)
        assert L == kernel

def test_cokernel_with_generators_and_projection():
    invariants, generators, projection = cokernel_with_generators_and_projection(
        3, [Vector([1,0,0]), Vector([0,0,1])])
    assert invariants == [0]
    assert generators == [Vector([0, 1, 0])]
    assert projection(Vector([123, 456, 789])) == [456]

    invariants, generators, projection = cokernel_with_generators_and_projection(
        3, [Vector([1,0,0]), Vector([0, 100, 0]), Vector([0,0,1])])
    assert invariants == [100]
    assert generators == [Vector([0, 1, 0])]
    assert projection(Vector([123, 456, 789])) == [56]

    invariants, generators, projection = cokernel_with_generators_and_projection(
        3, [Vector([0,2,2]), Vector([0,0,4])])
    assert invariants == [2, 4, 0]
    assert generators == [Vector([0, 1, 1]), Vector([0, 0, 1]), Vector([1, 0, 0])]
    assert projection(Vector([2, 3, 5])) == [1, 2, 2]
    assert projection(Vector([7, 7, 7])) == [1, 0, 7]

def test_cokernel_with_generators_and_projection_random():
    for N, R, _ in itertools.product(range(5), repeat=3):
        data = [Vector([random.randint(-10, 10) for _ in range(N)]) for _ in range(R)]
        invariants, generators, projection = cokernel_with_generators_and_projection(N, data)
        L = Lattice(N, data)
        L0 = L.copy()
        assert invariants == [d for d in L.nonzero_invariants() if d != 1] + [0] * (N - L.rank)
        for d, v in zip(invariants, generators):
            for k in range(1, d):
                assert k*v not in L
                assert k*v not in L0
            assert d*v in L
            assert d*v in L0
            L.add_vector(v)
        assert L.is_full()
        for i, gen in enumerate(generators):
            assert projection(gen) == [0] * i + [1] + [0] * (len(generators)-i-1), data
        zero = Vector.zero(N)
        for v in data:
            v1 = sum((x*gen for x, gen in zip(projection(v), generators)), start=zero)
            assert v1 - v in L0

def test_homology_with_generators_and_projection():
    invariants, generators, projection = homology_with_generators_and_projection(
        [Vector([0,0])],
        [Vector([7,7,7,7]),Vector([7,7,7,7])]
    )
    assert invariants == [0]
    assert generators == [Vector([1, -1])]
    assert projection(Vector([2,-2])) == [2]

    invariants, generators, projection = homology_with_generators_and_projection(
        [Vector([3,-2,0])],
        [Vector([2,2,2,2]),Vector([3,3,3,3]),Vector([0,0,0,0])]
    )
    assert invariants == [0]
    assert generators == [Vector([0, 0, 1])]
    assert projection(Vector([0, 0, 1])) == [1]
    assert projection(Vector([3, -2, 1])) == [1]
    assert projection(Vector([300, -200, 17])) == [17]


def test_homology_with_generators_and_projection_random():
    for K, N, R, _ in itertools.product(range(5), repeat=4):
        outgoing = [Vector([random.randint(-10,10) for _ in range(K)])
                    for _ in range(N)]
        kernel = relations_among(outgoing)
        incoming = [kernel.linear_combination(
                        Vector([random.randint(-10,10)
                                for _ in range(kernel.rank)]))
                    for _ in range(R)]
        invariants, generators, projection = homology_with_generators_and_projection(incoming, outgoing)
        L = Lattice(N, incoming)
        L0 = L.copy()
        assert invariants == [d for d in L.nonzero_invariants() if d != 1] + [0] * (kernel.rank - L.rank)
        for d, v in zip(invariants, generators):
            for k in range(1, d):
                assert k*v not in L
                assert k*v not in L0
            assert d*v in L
            assert d*v in L0
            L.add_vector(v)
        assert L == kernel
        for i, gen in enumerate(generators):
            assert projection(gen) == [0] * i + [1] + [0] * (len(generators)-i-1)
        zero = Vector.zero(N)
        for v in incoming:
            v1 = sum((x*gen for x, gen in zip(projection(v), generators)), start=zero)
            assert v1 - v in L0
