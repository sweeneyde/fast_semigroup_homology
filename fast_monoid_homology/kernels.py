"""
All of the functions here take in a list of mutable_lattice.Vector,
and return a basis for the set of relations between them.

If you give R vectors of length N, the result will be
a list of vectors of length R, the length of which is the nullity.
"""

from mutable_lattice import Lattice, Vector, relations_among, transpose

def mutable_lattice_kernel(vectors, *, verbose=False):
    if not vectors:
        return []
    if verbose:
        print(f"computing kernel of (R={len(vectors)})x(N={len(vectors[0])})")
    relations = relations_among(vectors).get_basis()
    if verbose:
        print(f"found {len(relations)} relations")
    return relations

def shuffled_columns(vectors):
    """Reshuffle all vectors in the same way, so that the least-used entries come first."""
    if not vectors:
        return
    N = len(vectors[0])
    scores = [0] * N
    rN = range(N)
    for v in vectors:
        for i in filter(v.__getitem__, rN):
            scores[i] += 10_000_000 + abs(v[i]).bit_length()
    columns_by_increasing_count = sorted(rN, key=scores.__getitem__)
    sort_action = [None] * N
    for i, col in enumerate(columns_by_increasing_count):
        sort_action[col] = i
    sort_action = Vector(sort_action)
    return [v.shuffled_by_action(sort_action) for v in vectors]

def vector_sortkey(vec):
    """A pair (-i, abs(vec[i])) indicating the first nonzero entry of the vector.
    Sorting vectors by this key attempts to make it as easy as possible to
    add these vectors to a given Lattice.
    """
    for i in filter(vec.__getitem__, range(len(vec))):
        return (-i, abs(vec[i]))
    return (-len(vec), 1)

def sorted_rows(vectors):
    return sorted(vectors, key=vector_sortkey)

def mutable_lattice_kernel_with_reshuffling(vectors, verbose=False):
    if not vectors:
        return []
    N = len(vectors[0])
    R = len(vectors)
    vectors = shuffled_columns(vectors)
    augmented_vectors = [
        Vector(v.tolist() + [0]*i + [1] + [0]*(R-1-i))
        for i, v in enumerate(vectors)
    ]
    augmented_vectors = sorted_rows(augmented_vectors)
    L = Lattice(N+R, augmented_vectors, maxrank=R)
    r2p = L._get_row_to_pivot()
    basis = L.get_basis()
    assert len(r2p) == len(basis)
    return [Vector(basis[i].tolist()[N:])
            for i in range(len(basis)) if r2p[i] >= N]

def shorten_vectors_to_hnf(vectors):
    N = len(vectors[0])
    R = len(vectors)
    columns = transpose(N, vectors)
    columns = sorted_rows(columns)
    columns = Lattice(R, columns, maxrank=len(columns)).get_basis()
    return transpose(R, columns)

def mutable_lattice_kernel_with_pre_hnf(vectors, *, verbose=False):
    if not vectors:
        return []
    if verbose:
        print("Shortening vectors with column ops...")
    vectors = shorten_vectors_to_hnf(vectors)
    return mutable_lattice_kernel_with_reshuffling(vectors, verbose=verbose)

def sage_kernel(vectors, *, verbose=False):
    # Only available if running the program from Sage instead of Python.
    if not vectors:
        return []
    N = len(vectors[0])
    R = len(vectors)
    from sage.all import Matrix, ZZ
    M = Matrix(ZZ, N, R, transpose(N, vectors))
    if verbose:
        from sage.misc.verbose import set_verbose
        set_verbose(2)
    K = M.right_kernel_matrix(
        # algorithm="default",
        algorithm="padic",
        # algorithm="pari",
        # basis="computed",
        # basis="echelon",
        basis="LLL",
    )
    if verbose:
        set_verbose(0)
    return [Vector(list(map(int, row))) for row in K]



