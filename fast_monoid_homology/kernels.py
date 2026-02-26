"""
All of the functions here take in a list of mutable_lattice.Vector,
and return a basis for the set of relations between them.

If you give R vectors of length N, the result will be
a list of vectors of length R, the length of which is the nullity.
"""

from mutable_lattice import relations_among, transpose, Vector, Lattice

def mutable_lattice_kernel(vectors, *, verbose=False):
    if not vectors:
        return []
    R = len(vectors)
    N = len(vectors[0])
    verbose = verbose or R > 1000
    if verbose:
        print(f"computing kernel of ({R=})x({N=})")
    relations = relations_among(vectors).get_basis()
    if verbose:
        print(f"{R}x{N} kernel found {len(relations)} relations")
    return relations

default_kernel = mutable_lattice_kernel

#######################################################
# The rest of this file has commented-out
# alternate implementations of kernels for comparison.
#######################################################

# def mutable_lattice_kernel_with_col_ops(vectors, *, verbose=False):
#     if not vectors:
#         return []
#     R = len(vectors)
#     N = len(vectors[0])
#     verbose = verbose or R > 1000
#     if verbose:
#         print(f"preworking columns for ({R=})x({N=}) kernel problem")
#     vectors = transpose(R, Lattice(R, transpose(N, vectors)).get_basis())
#     assert len(vectors) == R
#     if verbose:
#         print(f"shortened: ({N=}) --> (N={len(vectors[0])})")
#     relations = relations_among(vectors).get_basis()
#     if verbose:
#         print(f"{R}x{N} kernel found {len(relations)} relations")
#     return relations

# from cypari2 import Pari
# PARI = Pari()
# def _pari_hnf_kernel_with_flag(vectors, flag, verbose):
#     if not vectors:
#         return []
#     R = len(vectors)
#     N = len(vectors[0])
#     if verbose:
#         print(f"computing kernel of ({R=})x({N=})")
#     flat = []
#     for t in transpose(N, vectors):
#         flat.extend(t.tolist())
#     M = PARI.matrix(N, R, flat)
#     H, U, *maybe_P = M.mathnf(flag)
#     nullity = len(M) - len(H)
#     K = U[:nullity]
#     return [Vector(list(map(int, v))) for v in K]

# def pari_hnf_5_kernel(vectors, *, verbose=False):
#     return _pari_hnf_kernel_with_flag(vectors, 5, verbose)

# # from sage.all import Matrix as sage_Matrix, ZZ as sage_ZZ

# def sage_kernel_padic(vectors, *, verbose=False):
#     # Only available if running the program from Sage instead of Python.
#     if not vectors:
#         return []
#     if verbose:
#         print(f"computing kernel of (R={len(vectors)})x(N={len(vectors[0])})")
#     N = len(vectors[0])
#     R = len(vectors)
#     M = sage_Matrix(sage_ZZ, N, R, transpose(N, vectors))
#     K = M.right_kernel_matrix(
#         algorithm="default",
#         # algorithm="padic",
#         # algorithm="pari",
#         basis="computed",
#         # basis="echelon",
#         # basis="LLL",
#     )
#     if verbose:
#         print(f"found {K.ncols()} relations")
#     return [Vector(list(map(int, row))) for row in K]