"""
Here we use several techniques together to compute
cup products in the cohomology ring of a finite monoid.

The idea is to sandwich the chain complex B for B(M)
between two smaller chain complexes:

    Z(x)P ---f---> B ---p---> Q.

Here P is a hopefully-small ProjectiveResolution, tensored with Z,
and B is isomrphic to Z(x)E, where E is the bar resolution,
and we have a map P --> E lifting the identity, so that
by uniqueness of resolutions up to chain homotopy,
the first map f induces isomorphisms on (co)homology.

Then Q is the cellular chain complex for the quotient complex of B(M)
by under the collapsing scheme that Brown (1992) defines for
a monoid presented by a complete rewriting system,
so this quotient map p also induces isomorphisms on (co)homology.

Dualizing, we have the following cochain complexes with quasi-isos between:

    (Z(x)P)* <---f*--- B* <---p*--- Q*

Suppose we wish to compute the cup products from H^a x H^b --> H^(a + b).
* Generators for H^a(B) are given by p*(generators for H^a(Q*))
* Generators for H^b(B) are given by p*(generators for H^b(Q*))
* If cochains are represented in the cellular dual basis for B,
    the cup product is given by simply concatenating the tuples.
* We can compute generators for H^(a+b)(Z(x)P), and
    given a cup product cochain in B*, we can apply f*
    and project the result onto the given generators.
This lets us produce a 3D array:
    index 0: [generators alpha for H^a]
    index 1: [generators beta for H^b]
    index 2: [generators gamma for H^(a+b)]
    value: an integer describing how many times gamma appears in (alpha cup beta).

One issue is that computing the dual map p*
("which cells collapse onto this essential one")
is very intensive because there are many more cells in BM than in Q,
and the recursive nature of the collapsing scheme makes it difficult
to compute p* without just computing the entire matrix for p.

Instead, let {e1, ..., ei, ..., en} be the standard Z-basis for Z(x)P,
with dual basis {e1*, ..., ei*, ..., en*}.
Elements f(ei) in B* are exactly the data of of the map f.
We compute the coefficient of ei* in f*(p*(alpha) cup p*(beta)) as
    (f*(p*(alpha) cup p*(beta)))(ei)
    = (p*(alpha) cup p*(beta))(f(ei))
    = (p*(alpha) eating the left half of the simplices in f(ei))
        x (p*(beta) eating the right half of the simplices in f(ei))
    = alpha(p(left halves of f(ei))) x beta(p(right halves of f(ei)))

    ^ This is nice because p does not need to be computed for every simplex,
      only those present in the f(ei).

After computing this for each ei, we get a vector
in (Z(x)P)*, which we can project onto generators for the cohomology.
"""

from collections import Counter
from mutable_lattice import Vector
from ..projective_resolution import ProjectiveResolution
from .resolution_with_bar_comparison import ResolutionWithBarComparison
from .rewriting.crs import CRS
from .rewriting.op_to_crs import find_best_gens_crs

class FastCupComplex:
    """
    Store the data to compute cup products in the cohomology of a finite monoid
    """
    __slots__ = ["op", "res", "crs", "crs_elements",
                 "_crs_cohomology_generators_cache",
                 "_res_mapping_and_cohomology_cache",
                 "_cohomology_list_cache"]

    def __init__(self, argument, **find_best_gens_crs_kwargs):
        if isinstance(argument, list):
            op = argument
            for row in op:
                for x in row:
                    if not isinstance(x, int):
                        raise TypeError("multiplication table must be list[list[int]]")
            res = ResolutionWithBarComparison(op)
            crs = crs_elements = None
        elif isinstance(argument, ProjectiveResolution):
            op = argument.op
            res = ResolutionWithBarComparison(op, base_res=argument)
            crs = crs_elements = None
        elif isinstance(argument, ResolutionWithBarComparison):
            res = argument
            op = res.base_res.op
            crs = crs_elements = None
        elif isinstance(argument, CRS):
            crs = argument
            crs_elements = list(crs.elements())
            op = crs.op_table(crs_elements)
            res = ResolutionWithBarComparison(op)
        else:
            raise TypeError(f"Got unexpected type {type(argument)}")
        if crs is None:
            crs_elements, crs = find_best_gens_crs(op, **find_best_gens_crs_kwargs)
        self.op = op
        self.res = res
        self.crs = crs
        self.crs_elements = crs_elements
        self._crs_cohomology_generators_cache = {}
        self._res_mapping_and_cohomology_cache = {}
        self._cohomology_list_cache = []

    def cohomology_list(self, maxdim):
        if len(self._cohomology_list_cache) <= maxdim:
            cl = self.res.base_res.cohomology_list(maxdim)
            self._cohomology_list_cache = [
                list(Counter(cl_i).elements())[::-1]
                for cl_i in cl
            ]
        return self._cohomology_list_cache[:maxdim + 1]

    def cohomology(self, dim) -> dict[int, int]:
        if dim < len(self._cohomology_list_cache):
            return self._cohomology_list_cache[dim]
        else:
            return self.cohomology_list(dim)[dim]

    def crs_cohomology_with_generators(self, dim):
        cache = self._crs_cohomology_generators_cache
        if dim not in cache:
            cache[dim] = self.crs.cohomology_with_generators(dim)
        return cache[dim]

    def res_mapping_and_cohomology(self, dim):
        cache = self._res_mapping_and_cohomology_cache
        if dim not in cache:
            cache[dim] = self.res.mapping_and_cohomology(dim)
        return cache[dim]

    def cup_product_matrix(self, a, b):
        ca = self.cohomology(a)
        if not ca:
            return []
        cb = self.cohomology(b)
        if not cb:
            return [[] for _ in ca]
        cab = self.cohomology(a + b)
        if not cab:
            return [[[] for _ in cb] for _ in ca]
        a_invariants, a_gens = self.crs_cohomology_with_generators(a)
        assert a_invariants == ca
        b_invariants, b_gens = self.crs_cohomology_with_generators(b)
        assert b_invariants == cb
        e_images_in_bar, ab_cohomology = self.res_mapping_and_cohomology(a + b)
        elements = self.crs_elements
        cups_in_res = [[[0] * len(e_images_in_bar) for _ in b_gens] for _ in a_gens]
        for i_ab, e_image in enumerate(e_images_in_bar):
            for ab_tup_int, ab_tup_val in e_image.items():
                # Note: doing it in this order reduces
                # the number of calls to essential_representation.
                ab_tup = tuple(map(elements.__getitem__, ab_tup_int))
                a_tup_rep = self.crs.essential_representation(ab_tup[:a])
                b_tup_rep = self.crs.essential_representation(ab_tup[a:])
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
        assert ab_cohomology.get_invariants() == cab
        return matrix
