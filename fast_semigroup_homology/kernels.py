"""
All of the functions here take in a list of mutable_lattice.Vector,
and return a basis for the set of relations between them.

If you give R vectors of length N, the result will be
a list of vectors of length R, the length of which is the nullity.
"""

from mutable_lattice import (
    Vector, Lattice,
    relations_among,
    decompose_relations_among,
)

def mutable_lattice_c_verbose(vecs, /):
    if not vecs:
        return []
    R = len(vecs)
    N = len(vecs[0])
    keep = Vector(list(range(N)))
    vecs = [vec.shuffled_by_action(keep, N+R) for vec in vecs]
    for i, vec in enumerate(vecs):
        vec[N + i] = 1
    L = Lattice(N+R, maxrank=R)
    from tqdm import tqdm
    for vec in tqdm(vecs[::-1], miniters=1, ascii=True):
        L.add_vector(vec)
    basis = [Vector(v.tolist()[N:]) for v in L.get_basis()
             if next(filter(v.__getitem__, range(N + R))) >= N]
    return Lattice(R, basis, maxrank=len(basis)).get_basis()

def mutable_lattice_kernel(vectors, *, verbose=False):
    if not verbose:
        return relations_among(vectors).get_basis()
    if not vectors:
        return []
    R = len(vectors)
    N = len(vectors[0])
    print(f"getting kernel of ({R=})x({N=})...")
    relations, subproblem_rows, subproblems = decompose_relations_among(vectors)
    if verbose:
        print(f"found {len(relations)} easy relations")
        print(f"subproblems:", ",".join([f"(R={len(prob)})x(N={len(prob[0])})" for prob in subproblems]))
    for rows, subproblem in zip(subproblem_rows, subproblems, strict=True):
        if verbose:
            print(f"working on subproblem (R={len(subproblem)})x(N={len(subproblem[0])})...")
        subker = mutable_lattice_c_verbose(subproblem)
        if verbose:
            print(f"(R={len(subproblem)})x(N={len(subproblem[0])}) has {len(subker)} relations")
        for rel in subker:
            relations.append(rel.shuffled_by_action(rows, R))
    return Lattice(R, relations, maxrank=len(relations)).get_basis()

default_kernel = mutable_lattice_kernel

#######################################################
# The rest of this file has commented-out
# alternate implementations of kernels for comparison.
#######################################################

# from cypari2 import Pari
# PARI = Pari(4*1024*1024*1024)

# def _pari_hnf_kernel_with_flag(vectors, flag, verbose):
#     if not vectors:
#         return []
#     R = len(vectors)
#     N = len(vectors[0])
#     verbose = verbose or R > 1000
#     if verbose:
#         print(f"computing PARI kernel of ({R=})x({N=})")
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