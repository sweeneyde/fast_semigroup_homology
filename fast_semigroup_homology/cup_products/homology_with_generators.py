from mutable_lattice import Vector, Lattice, transpose, xgcd, relations_among

__all__ = ['homology_with_generators']

def _hnf(vectors):
    """Apply row operations to convert to HNF"""
    if not vectors:
        return vectors
    return Lattice(len(vectors[0]), vectors, maxrank=len(vectors)).get_basis()

def _augment_id(vectors):
    """Given a matrix A of row vectors, get the HNF of the augmented matix [A | id]"""
    R = len(vectors)
    return [Vector(v.tolist() + [0]*i + [1] + [0]*(R-i-1))
            for i, v in enumerate(vectors)]

def _inverse(vectors):
    """Get the inverse of an integer matrix."""
    N = len(vectors)
    augmented = _hnf(_augment_id(vectors))
    for i, v in enumerate(augmented):
        for j in range(N):
            if v[j] and i != j:
                raise ValueError("Not invertible")
    return [Vector(aug.tolist()[N:]) for aug in augmented]

def _smith_with_coefficients_unnormalized(N0, vectors):
    """Find a diagonal form of a matrix, along with the matrix of row operations used.

    Given a matrix A of row vectors, apply row operations to [A | id],
    along with column operations to A only, to get a matrix [D | S],
    where D=S*A*T is diagonal (but not necessarily square),
    and T (not provided) is the matrix of column operations used.
    Return the diagonal entries of D, along with the row vectors of S.
    """
    R = len(vectors)
    # Try column ops first because they're easier
    columns = transpose(N0, vectors)
    columns = _hnf(columns)
    N = len(columns)
    del vectors
    augmented = _hnf(_augment_id(transpose(R, columns)))
    # Repeat column ops and row ups until diagonal. Because we did the initial column steps,
    # the first N columns already have rank N, so we only need to check to the right
    # of the main diagonal.
    while any(any(map(aug.__getitem__, range(i + 1, N))) for i, aug in enumerate(augmented)):
        aug_columns = transpose(N+R, augmented)
        aug_columns[:N] = _hnf(aug_columns[:N])
        augmented = transpose(R, aug_columns)
        augmented = _hnf(augmented)
    assert N <= R
    invariants = [augmented[i][i] for i in range(N)] + [0] * (R - N)
    coefficient_vectors = [Vector(v.tolist()[N:]) for v in augmented]
    return invariants, coefficient_vectors

def _normalize_pair(d1, v1, d2, v2):
    """
    If [[d1,0],[0,d2]] is in SNF, with associated coefficient vectors v1 and v2,
    Fix the divisibility to produce the SNF [[gcd(d1,d2),0],[0,lcm(d1,d2)]],
    with the new associated coefficient vectors.
    """
    assert d1 and d2
    assert d2 % d1 != 0
    x, y, g = xgcd(d1, d2)
    # x and y satisfy:
    #   x * d1 + y * d2 == g
    #   (x//g)*d1 + (y//g)*d2 == 1
    # Implicitly do a column op, then apply an invertible matrix of row ops:
    # [    x      y] [d1  0|---v1---]  == [ g      d2*y |---------x*v1+y*v2----------]
    # [-d2//g d1//g] [d2 d2|---v2---]     [ 0  d1*d2//g |---(-d2//g)*d1+(d1//g)*v2---]
    # Since g divides d2 divides d2*y, an implicit column op restores diagonality.
    return g, x*v1+y*v2, d1//g*d2, (-d2//g)*v1+(d1//g)*v2

def _normalize_smith_with_coefficients(invariants, coefficient_vectors):
    """
    Repeatedly apply _normalize_pair to apply the SNF divisibility constraint
    to diagonal matrix with associated coefficient vectors
    """
    vecs, invars, zeros = [], [], []
    for d, v in zip(invariants, coefficient_vectors, strict=True):
        if d:
            assert d > 0
            invars.append(d)
            vecs.append(v)
        else:
            zeros.append(v)
    while True:
        perm = sorted(range(len(invars)), key=invars.__getitem__)
        invars = list(map(invars.__getitem__, perm))
        vecs = list(map(vecs.__getitem__, perm))
        changed = False
        for i in range(1, len(vecs)):
            d1, v1 = invars[i-1], vecs[i-1]
            d2, v2 = invars[i], vecs[i]
            if d2 % d1 != 0:
                d1, v1, d2, v2 = _normalize_pair(d1, v1, d2, v2)
                invars[i-1], vecs[i-1] = d1, v1
                invars[i], vecs[i] = d2, v2
                changed = True
        if not changed:
            return invars + [0] * len(zeros), vecs + zeros

def _smith_with_coefficients(N, vectors):
    """Find the SNF of a matrix, along with the matrix of row operations used.

    Given a matrix A of row vectors, apply row operations to [A | id],
    along with column operations to A only, to get a matrix [D | S],
    where D=S*A*T is in SNF, T (not provided) is the matrix of column operations used.
    Return the diagonal entries of D, along with the row vectors of S.
    """
    invariants, coefficient_vectors = _smith_with_coefficients_unnormalized(N, vectors)
    return _normalize_smith_with_coefficients(invariants, coefficient_vectors)

def _cokernel_with_generators_unfiltered(N, vectors):
    """
    Given vectors in Z^N, produce Z^N/(Z-span of vectors) as
    a list of cyclic summands (invariant factors), along with a vector in Z^N
    representing each cyclic summand. Include trivial summands Z/1.
    """
    vectors = _hnf(vectors)
    R = len(vectors)
    columns = transpose(N, vectors)
    invariants, column_coefficient_vectors = _smith_with_coefficients(R, columns)
    return invariants, _inverse(transpose(N, column_coefficient_vectors))

def cokernel_with_generators(N, vectors):
    """
    Given vectors in Z^N, produce Z^N/(Z-span of vectors) as
    a list of cyclic summands (invariant factors), along with a vector in Z^N
    representing each cyclic summand. Do NOT include trivial summands Z/1.
    """
    invariants, generators = _cokernel_with_generators_unfiltered(N, vectors)
    ones = sum(1 for d in invariants if d == 1)
    assert invariants[:ones] == [1] * ones
    return invariants[ones:], generators[ones:]

def homology_with_generators(incoming, outgoing):
    """Compute the homology at Z^R ---> Z^N ---> Z^K.

    `incoming` is a list of length R specifying the image in Z^N of each basis element of Z^R.
    `outgoing` is likewise a list of length N of vectors of length K.
    """
    kernel = relations_among(outgoing)
    relative = list(map(kernel.coefficients_of, incoming))
    invariants, rel_generators = cokernel_with_generators(kernel.rank, relative)
    generators = list(map(kernel.linear_combination, rel_generators))
    return invariants, generators
