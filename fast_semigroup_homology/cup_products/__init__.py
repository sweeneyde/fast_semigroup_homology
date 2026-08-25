from mutable_lattice import Vector, Lattice, transpose
from itertools import product
from math import gcd

from .resolution_with_bar_comparison import ResolutionWithBarComparison
from .homology_with_generators import homology_with_generators

def base_to_number(tup, n):
    res = 0
    for x in tup:
        res = res * n + x
    return res

def number_to_base(index, n, length):
    result = []
    for _ in range(length):
        index, rem = divmod(index, n)
        result.append(rem)
    assert index == 0
    result.reverse()
    return result

class BarComplex:
    __slots__ = ["op",
                 "identity",
                 "non_identity_elements",
                 "index_among_non_identity",
                 "cached_coboundary_matrices",
                 "cached_cohomology"]

    def __init__(self, op : list[list[int]]):
        self.op = op
        n = len(op)
        rn = range(n)
        [identity] = [e for e in rn
                    if all(op[x][e] == x == op[e][x] for x in rn)]
        self.identity = identity
        self.non_identity_elements = list(rn[:identity]) + list(rn[identity + 1:])
        self.index_among_non_identity = list(rn[:identity]) + [None] + list(range(identity, n - 1))
        self.cached_coboundary_matrices = {-1 : []}
        self.cached_cohomology = {}

    def tuple_to_index(self, tup):
        digits = list(map(self.index_among_non_identity.__getitem__, tup))
        return base_to_number(digits, len(self.op) - 1)

    def index_to_tuple(self, index, dim):
        digits = number_to_base(index, len(self.op) - 1, dim)
        return tuple(map(self.non_identity_elements.__getitem__, digits))

    def boundary_matrix(self, dim):
        if dim == 0:
            return [Vector([])]
        op = self.op
        t2i = self.tuple_to_index
        identity = self.identity
        boundary = []
        for tup in product(self.non_identity_elements, repeat=dim):
            face_vector = Vector.zero((len(op) - 1)**(dim - 1))
            face_vector[t2i(tup[1:])] += (-1)**0
            for i in range(1, dim):
                p = op[tup[i-1]][tup[i]]
                if p == identity:
                    continue
                face_vector[t2i(tup[:i-1] + (p,) + tup[i+1:])] += (-1)**i
            face_vector[t2i(tup[:-1])] += (-1)**dim
            boundary.append(face_vector)
        return boundary

    def coboundary_matrix(self, dim):
        if dim in self.cached_coboundary_matrices:
            return self.cached_coboundary_matrices[dim]
        boundary = self.boundary_matrix(dim + 1)
        cobounary = transpose((len(self.op) - 1)**dim, boundary)
        self.cached_coboundary_matrices[dim] = cobounary
        return cobounary

    def cohomology_with_generators(self, dim):
        if dim in self.cached_cohomology:
            return self.cached_cohomology[dim]
        incoming = self.coboundary_matrix(dim - 1)
        outgoing = self.coboundary_matrix(dim)
        invariants, generators = homology_with_generators(incoming, outgoing)
        assert 1 not in invariants
        rndim = range((len(self.op) - 1)**dim)
        bar_gens = [{self.index_to_tuple(i, dim) : gen[i]
                     for i in filter(gen.__getitem__, rndim)}
                    for gen in generators]
        result = invariants, bar_gens
        self.cached_cohomology[dim] = result
        return result

    def cohomology_freepart_with_generators(self, dim):
        invariants, bar_gens = self.cohomology_with_generators(dim)
        return [gen for d, gen in zip(invariants, bar_gens) if d == 0]

    def cup_product_exponents(self, a, b):
        a_invariants, a_gens = self.cohomology_with_generators(a)
        b_invariants, b_gens = self.cohomology_with_generators(b)
        if not a_gens:
            return a_invariants, b_invariants, []
        if not b_gens:
            return a_invariants, b_invariants, [[]] * len(a_gens)
        cob = self.coboundary_matrix(a + b - 1)
        coboundaries = Lattice((len(self.op) - 1)**(a + b),
                               cob,
                               maxrank=len(cob) + 1)
        result = []
        for da, a_gen in zip(a_invariants, a_gens):
            result_row = []
            for db, b_gen in zip(b_invariants, b_gens):
                cup = Vector.zero((len(self.op) - 1)**(a + b))
                for tup_a, ka in a_gen.items():
                    for tup_b, kb in b_gen.items():
                        cup[self.tuple_to_index(tup_a + tup_b)] += ka*kb
                g = gcd(da, db)
                if g:
                    for exponent in range(1, g + 1):
                        if exponent * cup in coboundaries:
                            result_row.append(exponent)
                            break
                    else:
                        raise AssertionError
                else:
                    with_cup = coboundaries.copy()
                    with_cup.add_vector(cup)
                    if with_cup.rank > coboundaries.rank:
                        result_row.append(0)
                    else:
                        exponent = 1
                        while exponent * cup not in coboundaries:
                            exponent += 1
                        result_row.append(exponent)
            result.append(result_row)
        return a_invariants, b_invariants, result

def cup_product_matrix(op, a, b):
    # Recall the universal coefficients theorem:
    bar = BarComplex(op)
    a_gens = bar.cohomology_freepart_with_generators(a)
    b_gens = bar.cohomology_freepart_with_generators(b)
    del bar
    sum_homology = ResolutionWithBarComparison(op).homology_freepart_generators_in_bar(a + b)
    return [
        [
            [
                sum(coeff * a_gen.get(tup[:a], 0) * b_gen.get(tup[a:], 0)
                    for tup, coeff in ab_gen.items())
                for ab_gen in sum_homology
            ]
            for b_gen in b_gens
        ]
        for a_gen in a_gens
    ]

