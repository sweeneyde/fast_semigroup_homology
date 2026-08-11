from mutable_lattice import Vector, transpose
from itertools import product

from .resolution_with_bar_comparison import ResolutionWithBarComparison
from .homology_with_generators import homology_with_generators

__all__ = ["cup_product_matrix"]

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

def cohomology_freepart_generators_in_bar(op, dim):
    if dim == 0:
        return [{(): 1}]
    assert dim >= 0
    n = len(op)
    rn = range(n)
    [identity] = [e for e in rn
                  if all(op[x][e] == x == op[e][x] for x in rn)]
    non_identity_elements = list(rn[:identity]) + list(rn[identity + 1:])
    index_among_non_identity = list(rn[:identity]) + [None] + list(range(identity, n - 1))

    def index_to_tuple(index, d):
        digits = number_to_base(index, n - 1, d)
        return tuple(non_identity_elements[x] for x in digits)

    def tuple_to_index(tup):
        digits = [index_among_non_identity[y] for y in tup]
        return base_to_number(digits, n - 1)

    def make_matrix(d):
        for tup in product(non_identity_elements, repeat=d):
            face_vector = Vector.zero((n - 1)**(d - 1))
            face_vector[tuple_to_index(tup[1:])] += (-1)**0
            face_vector[tuple_to_index(tup[:-1])] += (-1)**d
            for i in range(1, d):
                p = op[tup[i-1]][tup[i]]
                if p == identity:
                    continue
                face_vector[tuple_to_index(tup[:i-1] + (p,) + tup[i+1:])] += (-1)**i
            yield face_vector

    outgoing = list(make_matrix(dim))
    incoming = list(make_matrix(dim + 1))
    invariants, generators = homology_with_generators(
        transpose((n - 1)**(dim - 1), outgoing),
        transpose((n - 1)**dim, incoming)
    )
    result = []
    rndim = range((n - 1)**dim)
    for d, gen in zip(invariants, generators):
        if d == 0:
            result.append({
                index_to_tuple(i, dim) : gen[i]
                for i in filter(gen.__getitem__, rndim)
            })
    return result

def cup_product_matrix(op, a, b):
    cohomologies = {
        d: cohomology_freepart_generators_in_bar(op, d)
        for d in {a, b}
    }
    sum_homology = ResolutionWithBarComparison(op).homology_freepart_generators_in_bar(a + b)
    return [
        [
            [
                sum(coeff * a_gen.get(tup[:a], 0) * b_gen.get(tup[a:], 0)
                    for tup, coeff in ab_gen.items())
                for ab_gen in sum_homology
            ]
            for b_gen in cohomologies[b]
        ]
        for a_gen in cohomologies[a]
    ]

