"""
Use results related to Quillen's Theorem A
to prove that certain monoid maps f : M --> N
induce a homtopy equivalence Bf : BM --> BN
between classifying spaces.

The hypothesis of QTA is that the classifying space
of the undercategory/slice-category/comma-category
*/f is contractible (for each object *, but there's only one object * in M).
This undercategory has:
    objects are elements of n
    arrows n --> n' are elements m of M such that f(m)n=n'
So */f is like a "Cayley graph", showing how
elements of M translate via f between elements of N.

The classifying space of */f is thus the two-sided bar complex
B(*, M, N), with k-simplices written as [mk|...|m1]n0, so
    H_i(*/f)
    = H_i(chains in BM (tensor over ZM) ZN)
    = Tor^ZM_i(Z, ZN)
    = H_i(M; ZN).
which is exactly the sort of thing this software can compute.

We'll actually use a more refined statement from page 189
of Homological Algebra by Cartan and Eilenberg:
    f induces homology isomorphisms with all coefficients if:
        (1) f is a surjection
        (2) H_i(M; ZN) for all n > 0

It's well-known (e.g. Hatcher's Exercise 12 after Section 4.2)
that if we induce isomorphisms on homology with all coefficients
and on pi_1, then we induce a homotopy equivalence,
so detecting equivalences amounts checking
    (1) surjectivitity of f,
    (2) H_i(M; ZN)=0 for all n > 0
    (3) f induces isomorphisms on pi_1.

See also page 237 of:
    Marc Nunes. Cohomological results in monoid and category theory via classifying spaces.
    Journal of Pure and Applied Algebra 101 (1995) 213-244.

"""

from mutable_lattice import Vector
from . import ProjectiveResolution
from .group_completion import map_induces_isomorphism_on_group_completion

def assert_monoid_homomorphism(op1, op2, f):
    [e1] = [e for e in range(len(op1)) if all(op1[x][e] == x == op1[e][x] for x in range(len(op1)))]
    [e2] = [e for e in range(len(op2)) if all(op2[x][e] == x == op2[e][x] for x in range(len(op2)))]
    if f[e1] != e2:
        raise ValueError(f"f[identity={e1}] was {f[e1]}, expected {e2}")
    for x in range(len(op1)):
        for y in range(len(op1)):
            if f[op1[x][y]] != op2[f[x]][f[y]]:
                raise ValueError(f"f[op1[{x}][{y}]]=f[{op1[x][y]}]={f[op1[x][y]]}, but "
                                 f"op2[f[{x}]][f[{y}]]=op2[{f[x]}][{f[y]}]={op2[f[x]][f[y]]}")

def slice_resolution(op1, op2, f):
    assert_monoid_homomorphism(op1, op2, f)
    action = [Vector(op2[fx]) for fx in f]
    return ProjectiveResolution(op1, left_S_set_action=action)

def attempt_detect_equivalence(op1, op2, f):
    """
    Given a monoid homomorphism f : M --> N,
    determine if we can use Quillen's Theorem A to guarantee
    that f induces a homotopy equivalence between classifying spaces.
    """
    assert_monoid_homomorphism(op1, op2, f)
    return (
        set(f) == set(range(len(op2)))
        and map_induces_isomorphism_on_group_completion(op1, op2, f)
        and slice_resolution(op1, op2, f).has_trivial_homology()
    )
