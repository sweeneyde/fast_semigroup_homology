import itertools
from collections import Counter
from time import perf_counter
from datetime import timedelta, datetime
import multiprocessing
from hashlib import sha256
from pathlib import Path

try:
    from sage.all import (
        DeltaComplex,
        gap,
    )
except ImportError:
    print(f"Run this file with sage, not just Python:\n\n    sage {__file__}")
    quit(1)

from fast_semigroup_homology import ProjectiveResolution
from mutable_lattice import Vector, transpose
from cypari2 import Pari
PARI = Pari(4*1024*1024*1024)


if str(gap('LoadPackage("hap")')) != "true":
    print('Could not load GAP package "hap". Make sure it is installed.')
    quit(1)


RESULTS_FILE = Path(__file__).parent / "results.txt"

def fprint(*args):
    print(*args)
    with open(RESULTS_FILE, "a") as f:
        print(*args, file=f)

# These are extracted from the semisearch results folder
# monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units
MONOID_STRINGS = [
    "0",
    "01;10",
    "012;120;201",
    "0123;1230;2301;3012",
    "0123;1032;2301;3210",
    "01234;12340;23401;34012;40123",
    "01010;01011;23232;23233;01234",
    "012345;103254;254103;345012;430521;521430",
    "012345;103254;234501;325410;450123;541032",
    "010100;010101;232322;232323;010104;012345",
    "010100;010101;232322;232323;010144;012345",
    "010100;010111;232322;232333;012345;012354",
    "010100;010111;232322;232333;012345;230154",
    "010101;010110;232323;232332;012345;230154",
    "0123456;1234560;2345601;3456012;4560123;5601234;6012345",
    "0101000;0101001;2323222;2323223;0101004;0101005;0123456",
    "0101010;0101011;2323232;2323233;0101014;0101015;0123456",
    "0101010;0101011;2323232;2323233;0101014;2323235;0123456",
    "0101000;0101001;2323222;2323223;0101004;0101045;0123456",
    "0101000;0101001;2323222;2323223;0101004;0101055;0123456",
    "0101000;0101001;2323222;2323223;0101004;0101455;0123456",
    "0101000;0101001;2323222;2323223;0101004;2323255;0123456",
    "0101000;0101001;2323222;2323223;0101044;0101455;0123456",
    "0101000;0101001;2323222;2323223;0101044;0123055;0123456",
    "0101000;0101001;2323222;2323223;0101044;2323255;0123456",
    "0101010;0101011;2323232;2323233;0101014;2323255;0123456",
    "0101000;0101001;2323222;2323223;0101444;0101445;0123456",
    "0101000;0101001;2323222;2323223;0101454;0101545;0123456",
    "0101000;0101001;2323222;2323223;0101404;0101055;0123456",
    "0101000;0101001;2323222;2323223;0101404;2323255;0123456",
    "0101000;0101001;2323222;2323223;0101444;0101455;0123456",
    "0101000;0101001;2323222;2323223;0101444;0101555;0123456",
    "0101000;0101001;2323222;2323223;0101444;2323555;0123456",
    "0101000;0101001;2323222;2323223;0123444;0123555;0123456",
    "0101010;0101011;2323232;2323233;0101414;2323255;0123456",
    "0101000;0101011;2323222;2323233;0101044;0123456;0123465",
    "0101000;0101011;2323222;2323233;0101444;0123456;0123465",
    "0101000;0101111;2323222;2323333;0123456;0123564;0123645",
    "0120120;0120121;0120122;3453453;3453454;3453455;0123456",
]

MONOID_OPS = [
    [[int(x) for x in row]
     for row in s.split(";")]
    for s in MONOID_STRINGS
]

for op in MONOID_OPS:
    # Workaround for a bug in HAP: the identity must come first.
    # See https://github.com/gap-packages/hap/issues/129
    [e] = [e for e in range(len(op))
           if all(op[e][x] == x == op[x][e] for x in range(len(op)))]
    if e != 0:
        f = list(range(len(op)))
        f[0] = e
        f[e] = 0
        op[:] = [[f[op[fi][fj]] for fj in f] for fi in f]

##############################################
###   Alternate homology implementations   ###
##############################################

def sage_delta_complex_from_op(op, up_to_dimension):
    dimension_to_cells = {
        dim: tuple(itertools.product(range(len(op)), repeat=dim))
        for dim in range(0, up_to_dimension + 1)
    }
    cell_to_index = {
        cell: i
        for dim, cell_list in dimension_to_cells.items()
        for i, cell in enumerate(cell_list)
    }
    data = {0: [(),], 1: [(0, 0)]*len(op)}
    for dim in range(2, up_to_dimension + 1):
        cells = dimension_to_cells[dim]
        dimension_data = []
        for cell in cells:
            cell_data = []
            cell_data.append(cell_to_index[cell[1:]])
            for i in range(1, len(cell)):
                cell_data.append(cell_to_index[cell[:i-1] + (op[cell[i-1]][cell[i]],) + cell[i+1:]])
            cell_data.append(cell_to_index[cell[:-1]])
            dimension_data.append(tuple(cell_data))
        data[dim] = dimension_data
    return DeltaComplex(data=data, check_validity=False)

def sage_naive_homology(op, up_to_dimension):
    X = sage_delta_complex_from_op(op, up_to_dimension + 1)
    return [
        dict(Counter(map(int, reversed(X.homology(i).invariants()))))
        for i in range(1, up_to_dimension + 1)
    ]

def make_gap_bar_complex(op, up_to_dimension):
    one_indexed = [[x+1 for x in row] for row in op]
    gap_monoid = gap.MonoidByMultiplicationTable(one_indexed)
    return gap.BarComplexOfMonoid(gap_monoid, up_to_dimension)

def hap_naive_homology(op, up_to_dimension):
    gap_bar_complex = make_gap_bar_complex(op, up_to_dimension + 1)
    return [
        dict(Counter(reversed(
            list(map(int, gap.Homology(gap_bar_complex, i)))
        )))
        for i in range(1, up_to_dimension + 1)
    ]

def hap_contracted_homology(op, up_to_dimension):
    gap_bar_complex = make_gap_bar_complex(op, up_to_dimension + 1)
    gap_contracted_complex = gap.ContractedComplex(gap_bar_complex)
    return [
        dict(Counter(reversed(
            list(map(int, gap.Homology(gap_contracted_complex, i)))
        )))
        for i in range(1, up_to_dimension + 1)
    ]


############################################
###   Alternate kernel implementations   ###
############################################

def _pari_hnf_kernel_with_flag(vectors, flag):
    if not vectors:
        return []
    R = len(vectors)
    N = len(vectors[0])
    flat = []
    for t in transpose(N, vectors):
        flat.extend(t.tolist())
    M = PARI.matrix(N, R, flat)
    H, U, *maybe_P = M.mathnf(flag)
    nullity = len(M) - len(H)
    K = U[:nullity]
    return [Vector(list(map(int, v))) for v in K]

def pari_hnf_1_kernel(vectors, verbose=False):
    return _pari_hnf_kernel_with_flag(vectors, 1)

def pari_hnf_4_kernel(vectors, verbose=False):
    return _pari_hnf_kernel_with_flag(vectors, 4)

def pari_hnf_5_kernel(vectors, verbose=False):
    return _pari_hnf_kernel_with_flag(vectors, 5)

def pari_matkerint(vectors, verbose=False):
    if not vectors:
        return []
    R = len(vectors)
    N = len(vectors[0])
    flat = []
    for t in transpose(N, vectors):
        flat.extend(t.tolist())
    M = PARI.matrix(N, R, flat)
    K = M.matkerint()
    return [Vector(list(map(int, v))) for v in K]

from sage.all import Matrix as sage_Matrix, ZZ as sage_ZZ

def _sage_kernel_with_algorithm(vectors, algorithm):
    # Only available if running the program from Sage instead of Python.
    if not vectors:
        return []
    N = len(vectors[0])
    R = len(vectors)
    M = sage_Matrix(sage_ZZ, N, R, transpose(N, vectors))
    K = M.right_kernel_matrix(
        algorithm=algorithm,
        basis="computed",
    )
    return [Vector(list(map(int, row))) for row in K]

def sage_kernel_default(vectors, verbose=False):
    return _sage_kernel_with_algorithm(vectors, "default")

def sage_kernel_flint(vectors, verbose=False):
    return _sage_kernel_with_algorithm(vectors, "flint")

def sage_kernel_pari(vectors, verbose=False):
    return _sage_kernel_with_algorithm(vectors, "pari")

def sage_kernel_padic(vectors, verbose=False):
    return _sage_kernel_with_algorithm(vectors, "padic")

############################################################################
###   Homology functions that use the alternate kernel implementations   ###
############################################################################

def fast_homology(op, up_to_dimension):
    return ProjectiveResolution(op, check=False).homology_list(up_to_dimension)[1:]

def fast_homology_pari1(op, up_to_dimension):
    k = pari_hnf_1_kernel
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_pari4(op, up_to_dimension):
    k = pari_hnf_4_kernel
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_pari5(op, up_to_dimension):
    k = pari_hnf_5_kernel
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_parimatkerint(op, up_to_dimension):
    k = pari_matkerint
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_sage_default(op, up_to_dimension):
    k = sage_kernel_default
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_sage_flint(op, up_to_dimension):
    k = sage_kernel_flint
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_sage_pari(op, up_to_dimension):
    k = sage_kernel_pari
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

def fast_homology_sage_padic(op, up_to_dimension):
    k = sage_kernel_padic
    return ProjectiveResolution(op, check=False, kernel_implementation=k).homology_list(up_to_dimension)[1:]

#################################################
###   Main functions to do the benchmarking   ###
#################################################

def bench(homology_function, maxdim):
    total_time = 0
    num_iterations = 0
    while True:
        t0 = perf_counter()
        results = [homology_function(op, maxdim) for op in MONOID_OPS]
        t1 = perf_counter()
        total_time += t1 - t0
        num_iterations += 1
        if total_time > 10:
            break
    dt = timedelta(seconds=total_time/num_iterations)
    h = sha256(str(results).encode("ascii")).hexdigest()[:8]
    fprint(f"{maxdim=}, {dt}, ({h})")

def bench_with_timeout(homology_function, maxdim, timeout=10*60):
    p = multiprocessing.Process(target=bench, args=(homology_function, maxdim))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.terminate()
        fprint(f"{maxdim=} --> TIMEOUT after {timedelta(seconds=timeout)}")
        raise TimeoutError

HOMOLOGY_FUNCTIONS = [
    sage_naive_homology,
    hap_naive_homology,
    hap_contracted_homology,
    fast_homology,
    fast_homology_pari1,
    fast_homology_pari4,
    fast_homology_pari5,
    fast_homology_parimatkerint,
    fast_homology_sage_default,
    fast_homology_sage_flint,
    fast_homology_sage_padic,
    fast_homology_sage_pari,
]

def main():
    fprint(datetime.now())
    multiprocessing.set_start_method("spawn")

    fprint(f"===== BENCHMARKING {len(MONOID_OPS)} MONOIDS =====")

    for func in HOMOLOGY_FUNCTIONS:
        fprint(f"--- {func.__name__} ---")
        for maxdim in range(1, 20+1):
            try:
                bench_with_timeout(func, maxdim)
            except TimeoutError:
                break

if __name__ == "__main__":
    main()
