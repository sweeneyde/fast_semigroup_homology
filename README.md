# Fast Semigroup Homology

This is a package for computing the integral homology of finite semigroups.
It is very fast:

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

## Installation

`python -m pip install git+https://github.com/sweeneyde/fast_semigroup_homology@main`

## Usage

To compute the homology of an individual semigroup with the command line,
run `python -m fast_semigroup_homology -i "01;10" -d 10` as above, replacing
the string that follows `-i` with the multiplication table for some semigroup.
An integer following `-d` indicates the maximum dimension to compute to.

Alternatively, you can use this library from python code:
```python
>>> from fast_semigroup_homology import fast_integral_semigroup_homology as h
>>> h([[0,1],[1,0]], 10) # Homology of C2
[{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}]
```

If you have a folder of HDF5 files like those produced by [Semisearch](https://github.com/sweeneyde/semisearch),
you can instead run:

```bash
python -m fast_semigroup_homology -f /PATH/TO/FOLDER -o MAXORDER -d MAXDIM -c NUM_CORES
```
This will make a folder called `results` in the current working directory with .hdf5 files
containing the homology groups `[H_0, ... H_MAXDIM]` of the semigroups
from the files order1.hdf5 through orderMAXORDER.hdf5 in the folder passed after `-f`.
If NUM_CORES is specified, that many processes will compute in parallel.
Beware that some calculations can require a lot of memory, so making NUM_CORES
too high can make the program crash. If you are running out of
memory, you can also try [adding swap space](https://askubuntu.com/questions/178712/how-to-increase-swap-space).

My results folder is also stored in this repository, which you can browse and download.
This is not included when pip installs this package.