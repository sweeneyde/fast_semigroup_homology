from .resolution_with_bar_comparison import ResolutionWithBarComparison

from .bar_complex import BarComplex

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

