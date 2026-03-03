# Fast Semigroup Homology

This is a Python package for computing the integral homology of finite semigroups.
It is very fast in most cases:

```
$ python -m fast_semigroup_homology -i "010100;010100;232322;232322;010100;010100" -d 12
H_0: Z
H_1: trivial
H_2: Z
H_3: Z^2
H_4: Z^4
H_5: Z^8
H_6: Z^16
H_7: Z^32
H_8: Z^64
H_9: Z^128
H_10: Z^256
H_11: Z^512
H_12: Z^1024
Elapsed (wall) time: 0:00:00.001716
```

The above even quickly finishes if you pass `-d 10_000` to compute the homology
group in dimensions up to 10,000!

## Installation:

Install with pip: `pip install fast_semigroup_homology`

## Usage

### Homology of one semigroup from the command line:
To compute the homology of an individual semigroup with the command line,
run `python -m fast_semigroup_homology -i "01;10" -d 10` as above, replacing
the string that follows `-i` with the multiplication table for some semigroup.
Such a multiplication table must be entered as strings of digits separated by semicolons,
all wrapped in quotes. Letters (as in hexadecimal) are also supported
to allow semigroups of size at most 10+26=36 to be passed in quotes after `-i`.
The integer following `-d` indicates the maximum dimension to compute to.

### As a Python Function
With this package installed, you can use the
`fast_integral_semigroup_homology` function from python code:
```pycon
>>> from fast_semigroup_homology import fast_integral_semigroup_homology as h
>>> h([[0,1],[1,0]], 10) # Homology of C2
[{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}]
```

The result of `fast_integral_semigroup_homology(table, maxdim)`
is a list of length `maxdim+1` of dictionaries representing
the multisets of invariant factors of the homology groups in dimensions 0 up through `maxdim`, inclusive.
For example `{0: 1}` represents one copy of the infinite cyclic group Z,
while `{2: 1}` represents one copy of the cyclic group Z/2Z.

### Homology of HDF5 files of semigroup tables:

If you have a folder of [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)
files like those produced by [Semisearch](https://github.com/sweeneyde/semisearch),
you can run:

```bash
python -m fast_semigroup_homology -f /PATH/TO/FOLDER -o MAXORDER -d MAXDIM -c NUM_CORES
```
This will make a folder called `results` in the current working directory with .hdf5 files
containing the homology groups `[H_0, ... H_MAXDIM]` of the semigroups
from the files order1.hdf5 through orderMAXORDER.hdf5 in the folder passed after `-f`.
If NUM_CORES is specified, that many processes will compute in parallel.

[My results folder](https://github.com/sweeneyde/fast_semigroup_homology/tree/main/results)
is also stored in the git repository for this package, which you can browse and download.
This is not included when pip installs this package.

### Refining an existing HDF5 calculation

If you've already computed homology from a folder of HDF5 files as above,
you can attempt to extend these computations with the `-r` flag:

```bash
python -m fast_semigroup_homology -f /PATH/TO/FOLDER -r ./results/EXISTING.hdf5 -d MAXDIM -c NUM_CORES -k KERNEL_BOUND
```

This creates a new file "refined_EXISTING.hdf5"
in which calculations from "EXISTING.hdf5" are re-attempted,
organized into bins of semigroups whose homology computed so far has been identical.
If the calculations are found to be distinct, then
the semigroups are separated into bins accordingly and the homology of these bins
is extended in the same way recursively.

If any calculation within a bin attempts to compute the kernel of
some integer matrix larger than `KERNEL_BOUND`, then the calculation fails
and the bin is not refined any further. Therefore no
list of homology groups produced by this refinement
process is a prefix of any other.

The final format of `refined_EXISTING.hdf5`
is identical to the input `EXISTING.hdf5`,
and a similar `refined_EXISTING.md` is generated to summarize
the calculation. Question mark emojis ("❓")
in the markdown files indicate computations that were abandoned
for exceeding `KERNEL_BOUND`.

## How It Works

We compute semigroup homology by first reducing to monoid homology,
adjoining a unit as necessary. Then, according to Theorem 2 of
(Adams, W. and Rieffel, M.
 Adjoint functors and derived functors with an application
 to the cohomology of semigroups. *J. Algebra*, 7 (1967), 25-34),
we can get the same homology groups by reducing from a monoid *S*
to the monoid *eSe* for some idempotent *e* of *S*,
so long as *eSe* is a left ideal or right ideal
of *S*.

After reducing to smaller monoids as applicable,
the computation computes the homology *H_i(S)*
as a Tor group Tor^Z[*S*]_i(Z, Z), where Z[*S*] is the monoid
ring for *S*, and Z is a Z[*S*]-module with trivial *S*-action.
To do this, we need a projective resolution for Z as a Z[*S*]-module.

Forming such a projective resolution is an alternating process of taking kernels
and covering those kernels with maps from a new projective Z[*S*]-module.
The projective modules we use are direct sums of modules of the form Z[*S*]*e*,
where *e* is some idempotent of *S*. Kernels of integer matrices
are computed using the implementation in
./fast_semigroup_homology/kernels.py, and those kernels
are then covered by maps from a new module using the implementation
in ./fast_semigroup_homology/find_generating_subset.py.

We apply two major additional optimizations in ./fast_semigroup_homology/projective_resolution.py:

1. We use the `mutable_lattice.Lattice.decompose()`
method to identify when a kernel splits along the Z[*S*]*e* summands,
so we can cover the two summands of the kernel independently,
thus "forking" the projective resolution into multiple branches of a tree.

2. We cache and re-use parts of the resolution already computed,
so that infinite resolutions can be represented as finite graphs with cycles.
This is what allows large calculations like the 1000th homology group
of the semigroup from the top of this README.

## Performance Notes

* Computing homology is extremely fast for many semigroups,
but much slower for others. The time spent calculating
homology for a big list of semigroups might be dominated
by a tiny handful of the most difficult examples. In
these most difficult cases, the complexity of
the computation is exponential in `MAXDIM`.

* The `fast_integral_semigroup_homology` function
aims to provide shortcuts and sensible default thresholds to minimize
the computation time expended finding small
projective resolutions. You can get more fine-grained control by using
the `ProjectiveResolution` class directly.

* It appears that the majority of calculation time occurs
in computing the kernels of integer matrices.
The `mutable_lattice.relations_among`
function proved to be the fastest for the data provided by actual
semigroups, which appear to usually small integer entries including many zeros.
However, `ProjectiveResolution` allows dependency-injecting your
own kernel function. See ./fast_semigroup_homology/kernels.py
for some commented-out examples of alternate implementations.
This dependency-injection also allows bounding the size of kernel computation
that is allowed, as is done in ./fast_semigroup_homology/refine_hdf5.py.

* The most intensive calculations tend to take lots of memory,
and this is especially relevant if you are running multiple
processes in parallel with the `-c` flag.
[Adding swap space](https://askubuntu.com/questions/178712/how-to-increase-swap-space)
can help offload some of the required memory onto your hard drive
and extend the computations your computer can complete without crashing.