from .projective_resolution import ProjectiveResolution

__all__ = ["fast_integral_monoid_homology"]

def permute_op(op, g, flip):
    g_inv = [None] * len(op)
    for i in range(len(op)):
        g_inv[g[i]] = i
    if flip:
        return [[g_inv[op[gj][gi]] for gj in g] for gi in g]
    else:
        return [[g_inv[op[gi][gj]] for gj in g] for gi in g]

def _easy_attempt(
        op,
        maxdim,
        *,
        peek_dimension=4,
        min_size_to_try_harder=20,
        verbose=False,
):
    # first, just try computing with a single resolution:
    fast_res = ProjectiveResolution(op)
    fast_res.extend_to_dimension(
        peek_dimension,
        sloppy_last_cover=False,
        max_size_to_ensure_minimal=-1,
        verbose=verbose)
    if fast_res.total_cost(peek_dimension) < min_size_to_try_harder:
        return fast_res.homology_list(
            maxdim,
            max_size_for_extra_greedy=-1,
            max_size_to_ensure_minimal=-1,
            verbose=verbose)
    return None

def _hard_attempt(
        op,
        *,
        peek_dimension=4,
        verbose=False,        
):
    # when that's getting too big, try shuffling things around
    # until a better resolution is available:
    rn = list(range(len(op)))
    permutations = [rn[:], rn[::-1]]
    op_list = [
        permute_op(op, g, flip)
        for g in permutations
        for flip in (False, True)
    ]
    res_list = [ProjectiveResolution(op) for op in op_list]
    for res in res_list:
        res.extend_to_dimension(
            peek_dimension,
            sloppy_last_cover=False,
            max_size_for_extra_greedy=100,
            max_size_to_ensure_minimal=100,
            verbose=verbose,
        )
    res = min(res_list, key=lambda r: r.total_cost(peek_dimension))
    return res

def integral_monoid_homology(
        op,
        maxdim,
        *,
        peek_dimension=None,
        min_size_to_try_harder=20,
        verbose=False,
):
    """
    Compute some integral monoid homology via a projective resolution.

    >>> integral_monoid_homology([[0,1],
    ...                           [1,0]], 10)
    [{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}]
    >>> integral_monoid_homology([[0,1,0,1,0],
    ...                           [0,1,0,1,1],
    ...                           [2,3,2,3,2],
    ...                           [2,3,2,3,3],
    ...                           [0,1,2,3,4]], 10)
    [{0: 1}, {}, {0: 1}, {}, {}, {}, {}, {}, {}, {}, {}]    
    >>> integral_monoid_homology([[0,1,0,1,0,0],
    ...                           [0,1,0,1,0,1],
    ...                           [2,3,2,3,2,2],
    ...                           [2,3,2,3,2,3],
    ...                           [0,1,0,1,0,4],
    ...                           [0,1,2,3,4,5]], 10)
    [{0: 1}, {}, {0: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1}]

    Strategy: First try any projective resolution you can find.

    If we start doing that and realize it'll be too hard, instead try some variations,
    and pick the resolution that looks the smallest so far before continuing to extend and use it.
    """
    if peek_dimension is None:
        peek_dimension = max(0, min(4, maxdim - 2))
    easy_result = _easy_attempt(op, maxdim, peek_dimension=peek_dimension, min_size_to_try_harder=min_size_to_try_harder, verbose=verbose)
    if easy_result is not None:
        return easy_result
    res = _hard_attempt(op, peek_dimension=peek_dimension, verbose=verbose)
    return res.homology_list(
        maxdim,
        max_size_for_extra_greedy=100,
        max_size_to_ensure_minimal=100,
        verbose=verbose,
    )

def maybe_adjoin_1(op):
    rn = range(len(op))
    for e in rn:
        if all(op[e][x] == x == op[x][e] for x in rn):
            return op
    new_op = [list(row)+[i] for i, row in enumerate(op)]
    new_op.append(list(range(len(op) + 1)))
    return new_op

def equivalent_submonoid(op):
    # Search for an idempotent e with eSe=eS or eSe=Se.
    # Replace S with eSe and repeat until done.
    while True:
        rn = range(len(op))
        E = [e for e in rn if op[e][e] == e]
        good_subsets = []
        for e in E:
            eS = {op[e][x] for x in rn}
            Se = {op[x][e] for x in rn}
            eSe = eS & Se
            if eSe == eS or eSe == Se:
                good_subsets.append(eSe)
        T = min(good_subsets, key=len)
        if len(T) == len(op):
            return op
        T_index = {t: i for i, t in enumerate(T)}
        T_op = [[T_index[op[Ti][Tj]] for Tj in T] for Ti in T]
        op = T_op

def fast_integral_semigroup_homology(
        op,
        maxdim,
        *,
        verbose=False,
):
    """Get the integral homology of a semigroup.
    First apply reductions to smaller semigroups/monoids if applicable.

    >>> fast_integral_semigroup_homology([[x]*100 for x in range(100)], 10)
    [{0: 1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    >>> fast_integral_semigroup_homology([[(x+y)%2 for x in range(100)] for y in range(100)], 10)
    [{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}]
    """
    if not any(op[0]):
        # Extremely common for semigroups: having a (left) zero element.
        return [{0: 1}] + [{}] * maxdim
    op = equivalent_submonoid(maybe_adjoin_1(op))
    if len(op) == 1:
        return [{0: 1}] + [{}] * maxdim
    elif len(op) in (2, 3):
        return [{0: 1}] + [{len(op): 1} if i % 2 == 1 else {} for i in range(1, maxdim)]
    return integral_monoid_homology(op, maxdim, verbose=verbose)
