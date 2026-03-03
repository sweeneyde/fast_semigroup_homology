Results go in this folder.

The subfolder names follow the names of the folders output
by [Semisearch](https://github.com/sweeneyde/semisearch).

Each job here produces a binary .hdf5 file with the full homology results
alongside a more human-readable markdown file summarizing the computations.

*   The .hdf5 file will have a dataset named `homology_group_lists`
that stores ascii representations of homology results, and several
datasets named `order1`, `order2`, etc, storing indexes into the `homology_group_lists`.
Consider using [HDFView](https://www.hdfgroup.org/download-hdfview/)
to browse the structure of HDF5 files.

* The corresponding .md markdown file will contain a table that represents the
`homology_group_lists` dataset with unicode characters like "ℤ²×𝐶₂",
along with the number of semigroups found with each list homology groups.
The infinite cyclic group is denoted with ℤ, finite cyclic groups are denoted with
𝐶 with subscripts, multiplicities are denoted with superscripts, and the trivial group
is denoted by "·".
