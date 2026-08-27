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
    del vectors
    columns = _hnf(columns)
    N = len(columns)
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
    to a diagonal matrix with associated coefficient vectors
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
    where D=S*A*T is in SNF, and T (not provided) is the matrix of column operations used.
    Return the diagonal entries of D, along with the row vectors of S.
    """
    invariants, coefficient_vectors = _smith_with_coefficients_unnormalized(N, vectors)
    return _normalize_smith_with_coefficients(invariants, coefficient_vectors)

def _smith_with_right_coefficients(N, vectors):
    """Find the SNF of a matrix, along with the matrix of column operations used.

    Given a matrix A of row vectors, apply row operation matrix S
    and column operation matrix T such that D=S*A*T is in SNF.
    Return the diagonal entries of D, along with the row vectors of T.
    """
    columns = transpose(N, vectors)
    R = len(vectors)
    invariants, coeff = _smith_with_coefficients(R, columns)
    return invariants, transpose(len(coeff), coeff)

class Cokernel:
    """Represent a quotient group Z^N/vectors"""
    __slots__ = [
        # The ambient dimension our vectors live in
        "N",
        # The vectors we're quotienting by, as HNF
        "vectors",
        # The invariants for the Z/d summands, including 0s (Z summands),
        #     but not 1s (trivial summands)
        "_invariants",
        # The generators of Z^N that get mapped to to the generators
        #     of the cyclic summand Z/d of Z^N/vectors
        "_generators",
        # The transposed invertible matrix T such that S(vectors)T = D is in Smith normal form.
        "_smith_T",
        # The image of each standard basis vector in the nontrivial summands
        "_standard_basis_to_nontrivial",
    ]

    def __init__(self, N : int, vectors : list[Vector]):
        self.N = N
        self.vectors = _hnf(vectors)
        self._invariants = None
        self._smith_T = None
        self._generators = None
        self._standard_basis_to_nontrivial = None

    def get_invariants(self) -> list[int]:
        """
        Get the invariants for the Z/d summands in this quotient,
        including 0s (Z summands), but not 1s (trivial summands).
        """
        if self._invariants is None:
            L = Lattice(self.N, self.vectors, maxrank=len(self.vectors))
            self._invariants = [d for d in L.invariants() if d != 1]
        return self._invariants

    def get_smith_T(self) -> list[Vector]:
        """
        Get the invertible matrix T such that S(vectors)T = D is in SNF.
        """
        if self._smith_T is None:
            all_invariants, smith_T = _smith_with_right_coefficients(self.N, self.vectors)
            invariants = [d for d in all_invariants if d != 1]
            if self._invariants is not None:
                assert self._invariants == invariants
            self._invariants = invariants
            self._smith_T = smith_T
        return self._smith_T

    def get_generators(self) -> list[Vector]:
        if self._generators is None:
            V = _inverse(self.get_smith_T())
            num_ones = self.N - len(self._invariants)
            self._generators = V[num_ones:]
        return self._generators

    def get_standard_basis_to_nontrivial(self) -> list[Vector]:
        if self._standard_basis_to_nontrivial is None:
            T = self.get_smith_T()
            num_ones = self.N - len(self._invariants)
            T_slice = transpose(self.N, transpose(self.N, T)[num_ones:])
            self._standard_basis_to_nontrivial = T_slice
        return self._standard_basis_to_nontrivial

    def projection(self, v : Vector) -> list[int]:
        s2nontrivial = self.get_standard_basis_to_nontrivial()
        invariants = self.get_invariants()
        result = Vector.zero(len(invariants))
        for vi, gen_coeffs_i in zip(v, s2nontrivial, strict=True):
            result += vi * gen_coeffs_i
        return [x % d if d else x for x, d in zip(result, invariants, strict=True)]

class SubQuotient:
    __slots__ = ["L", "relative_coker", "_generators"]
    def __init__(self, L, vectors):
        self.L = L
        self.relative_coker = Cokernel(L.rank, list(map(L.coefficients_of, vectors)))
        self._generators = None

    def get_generators(self):
        if self._generators is None:
            relative_generators = self.relative_coker.get_generators()
            self._generators = list(map(self.L.linear_combination, relative_generators))
        return self._generators

    def get_invariants(self):
        return self.relative_coker.get_invariants()

    def projection(self, v):
        return self.relative_coker.projection(self.L.coefficients_of(v))

def homology(incoming, outgoing):
    """Compute the homology at Z^R ---> Z^N ---> Z^K.

    `incoming` is a list of length R specifying the image in Z^N of each basis element of Z^R.
    `outgoing` is likewise a list of length N of vectors of length K.
    """
    return SubQuotient(relations_among(outgoing), incoming)
