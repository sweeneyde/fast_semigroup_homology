from mutable_lattice import Lattice, Vector
from fast_monoid_homology.kernels import (
    mutable_lattice_kernel,
    mutable_lattice_kernel_with_reshuffling,
    mutable_lattice_kernel_with_pre_hnf,
    sage_kernel,
)

ker_funcs = [
    mutable_lattice_kernel,
    mutable_lattice_kernel_with_reshuffling,
    mutable_lattice_kernel_with_pre_hnf,
]

import pytest
import random

@pytest.mark.parametrize("ker_func", ker_funcs)
@pytest.mark.parametrize("vectors,expected", [
    ([], []),
    ([Vector([1])], []),
    ([Vector([2])], []),
    ([Vector([0])], [Vector([1])]),
    ([Vector([5]), Vector([0])], [Vector([0, 1])]),
    ([Vector([1]), Vector([-1])], [Vector([1, 1])]),
    ([Vector([1]), Vector([1])], [Vector([1, -1])]),
    ([Vector([1]), Vector([2])], [Vector([2, -1])]),
    ([Vector([2]), Vector([3])], [Vector([3, -2])]),
    ([Vector([0]), Vector([0])], [Vector([1, 0]), Vector([0, 1])]),
    ([Vector([1]), Vector([2]), Vector([3])], [Vector([1, 1, -1]), Vector([0, 3, -2])]),
    ([Vector([1, 2])], []),
    ([Vector([0, -1])], []),
    ([Vector([0, 0])], [Vector([1])]),
    ([Vector([1, 2]), Vector([3, 4])], []),
    ([Vector([1, 2]), Vector([-1, -2])], [Vector([1, 1])]),
    ([Vector([1, 2]), Vector([3, 4]), Vector([4, 6])], [Vector([1, 1, -1])]),
    ([Vector([100, 200]), Vector([3, 4]), Vector([4, 6])], [Vector([1, 100, -100])]),
])
def test_kernels(ker_func, vectors, expected):
    assert ker_func(vectors) == expected


@pytest.mark.parametrize("ker_func", ker_funcs)
@pytest.mark.parametrize("N", [0, 1, 2, 3])
@pytest.mark.parametrize("R", [0, 1, 2, 3])
def test_kernels_random(ker_func, N, R):
    pool = [-10*100, -3, -2, -1, 0, 1, 2, 3, 10*100]
    random.seed(0)
    for _ in range(10):
        data = [random.choices(pool, k=N) for _ in range(R)]
        vectors = list(map(Vector, data))
        ker_basis = ker_func(vectors)
        assert [v.tolist() for v in vectors] == data
        for relation in ker_basis:
            s = Vector([0] * N)
            for coeff, vec in zip(relation, vectors, strict=True):
                s += coeff * vec
            assert s.tolist() == [0] * N, (vectors, ker_basis)
        rank = Lattice(N, vectors).rank
        nullity = len(ker_basis)
        assert rank + nullity == R
        # kernels are always saturated:
        assert Lattice(R, ker_basis).nonzero_invariants() == [1] * nullity