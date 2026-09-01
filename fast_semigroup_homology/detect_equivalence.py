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

from itertools import combinations
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

def equivalence_to_quotient(op):
    for op2, f in proper_nontrivial_quotient_semigroups(op):
        if attempt_detect_equivalence(op, op2, f):
            return op2, f
    return None

def proper_nontrivial_quotient_semigroups(op):
    n = len(op)
    for cong in congruences(op):
        adjacencies = [[] for _ in range(len(op))]
        for x, y in cong:
            adjacencies[x].append(y)
            adjacencies[y].append(x)
        components = set(map(frozenset, adjacencies))
        components = sorted(map(sorted, components))
        if len(components) == 1 or len(components) == n:
            continue
        f = [None] * len(op)
        for fx, comp in enumerate(components):
            for x in comp:
                f[x] = fx
        op2 = [[f[op[comp1[0]][comp2[0]]] for comp2 in components] for comp1 in components]
        yield op2, f

def congruences(op):
    # Use Bernhard Ganter's NextClosure algorithm,
    # Alg 1.7 from from https://doi.org/10.1007/978-3-642-11928-6_22
    n = len(op)
    rn = range(n)
    potential_pairs = [(i, j) for i in rn for j in rn if i < j]
    # A stack of (element_added, closure) pairs
    history = [((-1, -1), {(x, x) for x in rn})]
    while True:
        yield history[-1][1]
        # Add in the least significant thing we can
        for p in reversed(potential_pairs):
            if p in history[-1][1]:
                continue
            # Undo previous additions less significant than p,
            # so only things p necessitates get added.
            while p < history[-1][0]:
                history.pop()

            A = history[-1][1]
            Ap = set()
            for q in new_implied_congruences(op, A, p):
                assert q not in A
                if q < p:
                    # If adding this p added some more significant q,
                    # Don't count this now; count it later when q
                    # has its turn to be added.
                    break
                Ap.add(q)
            else:
                Ap.update(A)
                # No break: everything has waited its turn.
                history.append((p, Ap))
                break # back to yield
        else:
            # No break: nothing yielded
            assert len(history[-1][1]) == n * (n + 1) // 2
            return

def whiskered_congruences(op, pair):
    a, b = pair
    for x in range(len(op)):
        yield tuple(sorted((op[a][x], op[b][x])))
        yield tuple(sorted((op[x][a], op[x][b])))

def get_components(n, pairs):
    adjacencies = [[] for _ in range(n)]
    for x, y in pairs:
        adjacencies[x].append(y)
        adjacencies[y].append(x)
    done = set()
    for start in range(n):
        if start in done:
            continue
        component = []
        stack = [start]
        while stack:
            x = stack.pop()
            if x not in done:
                done.add(x)
                component.append(x)
                stack.extend(adjacencies[x])
        component.sort()
        yield component

def new_implied_congruences(op, pairs, new_pair):
    pairs = pairs.copy()
    for pair in whiskered_congruences(op, new_pair):
        if pair not in pairs:
            yield pair
            pairs.add(pair)
    while True:
        new_from_transitive = set()
        for comp in get_components(len(op), pairs.copy()):
            for pair in combinations(comp, 2):
                if pair not in pairs:
                    yield pair
                    new_from_transitive.add(pair)
                    pairs.add(pair)
        if not new_from_transitive:
            break
        progress = False
        for pair1 in new_from_transitive:
            for pair2 in whiskered_congruences(op, pair1):
                if pair2 not in pairs:
                    yield pair2
                    pairs.add(pair2)
                    progress = True
        if not progress:
            break
