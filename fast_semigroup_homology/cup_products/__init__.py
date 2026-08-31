from .resolution_with_bar_comparison import ResolutionWithBarComparison

from .bar_complex import BarComplex
from .rewriting.op_to_crs import find_best_gens_crs

from mutable_lattice import Vector


def cup_product_freepart_matrix(op, a, b):
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

def cup_product_matrix(op, a, b, **kwargs):
    elements, crs = find_best_gens_crs(op, **kwargs)
    a_invariants, a_gens = crs.cohomology_with_generators(a)
    b_invariants, b_gens = crs.cohomology_with_generators(b)
    res = ResolutionWithBarComparison(op)
    e_images_in_bar, ab_cohomology = res.mapping_and_cohomology(a + b)
    del res
    cups_in_res = [[[0] * len(e_images_in_bar) for _ in b_gens] for _ in a_gens]
    for i_ab, e_image in enumerate(e_images_in_bar):
        for ab_tup_int, ab_tup_val in e_image.items():
            ab_tup = tuple(map(elements.__getitem__, ab_tup_int))
            a_tup_rep = crs.essential_representation(ab_tup[:a])
            b_tup_rep = crs.essential_representation(ab_tup[a:])
            a_tup_values_on_ab_tup = [
                sum(a_tup_rep[a_tup1] * cochain_a[a_tup1]
                    for a_tup1 in a_tup_rep.keys() & cochain_a.keys())
                for cochain_a in a_gens
            ]
            b_tup_values_on_ab_tup = [
                sum(b_tup_rep[b_tup1] * cochain_b[b_tup1]
                    for b_tup1 in b_tup_rep.keys() & cochain_b.keys())
                for cochain_b in b_gens
            ]
            for i_a, a_value in enumerate(a_tup_values_on_ab_tup):
                for i_b, b_value in enumerate(b_tup_values_on_ab_tup):
                    cups_in_res[i_a][i_b][i_ab] += a_value * b_value * ab_tup_val
    matrix = [
        [ab_cohomology.projection(Vector(cup)) for cup in cup_row]
        for cup_row in cups_in_res
    ]
    ab_invariants = ab_cohomology.get_invariants()
    return a_invariants, b_invariants, ab_invariants, matrix

