from pathlib import Path

def check_associative(op):
    for row in op:
        if len(row) != len(op):
            raise ValueError("semigroup operation table was not square")
    rn = range(len(op))
    for row in op:
        for entry in row:
            if entry not in rn:
                raise ValueError(f"{len(op)}x{len(op)} semigroup operation table cannot contain '{entry:x}'")
    for i in rn:
        for j in rn:
            ij = op[i][j]
            for k in rn:
                jk = op[j][k]
                if op[ij][k] != op[i][jk]:
                    raise ValueError(f"op[op[{i}][{j}]][{k}]=op[{ij}][{k}]={op[ij][k]}, but "
                                        f"op[{i}][op[{j}][{k}]]=op[{i}][{jk}]={op[i][jk]}")

def main_individual(opstring, maxdim, verbose):
    op = [[int(c, 36) for c in line] for line in opstring.strip().split(";")]
    check_associative(op)
    from .homology import fast_integral_semigroup_homology
    import time
    import datetime
    t0 = time.time()
    hlist = fast_integral_semigroup_homology(op, maxdim, verbose=verbose)
    for i, h in enumerate(hlist):
        if not h:
            group_name = "trivial"
        else:
            group_name_parts = []
            for divisor, count in h.items():
                summand = "Z" if divisor == 0 else f"C{divisor}"
                assert count >= 1
                if count > 1:
                    summand += f"^{count}"
                group_name_parts.append(summand)
            group_name = " x ".join(group_name_parts)
        print(f"H_{i}: {group_name}")
    t1 = time.time()
    dt = datetime.timedelta(seconds=t1-t0)
    print(f"Elapsed (wall) time: {dt}")

def main():
    import argparse
    parser = argparse.ArgumentParser(prog="fast_monoid_homology",
                                     description="compute homology of finite semigroups",
                                     usage="")
    parser.add_argument("-f", "--folder",
                        help="folder of hdf5 files full of semigroups",
                        type=Path, default=None)
    parser.add_argument("-i", "--individual",
                        help="a multiplication table of an individual semigroup",
                        type=str, default=None)
    parser.add_argument("-v", "--verbose",
                        help="flag to display more of the computation process",
                        action='store_true')
    parser.add_argument("-c", "--num_cores",
                        help="the number of cores used to use in parallel",
                        default=None, type=int)
    parser.add_argument("-d", "--maxdim",
                        help="the maximum dimension of homology groups to compute",
                        default=6, type=int)
    parser.add_argument("-o", "--maxorder",
                        help="the maximum order of semigroup for which to compute the homology",
                        type=int)
    parser.add_argument("--chunksize",
                        help="the number of semigroups in a batch to send to a multiprocessing worker",
                        type=int, default=1)
    args = parser.parse_args()

    if args.individual is not None:
        if args.num_cores is not None:
            raise ValueError(f"Can't pass -c {args.num_cores} with -i.")
        if args.folder is not None:
            raise ValueError(f"Can't pass -f {args.folder} with -i.")
        if args.maxorder is not None:
            raise ValueError(f"Can't pass -o {args.maxorder} with -i.")
        main_individual(
            opstring=args.individual,
            maxdim=args.maxdim,
            verbose=args.verbose
        )
    else:
        if args.num_cores is None:
            args.num_cores = 4
        if args.folder is None:
            raise ValueError("Must pass either -i or -f")
        if args.maxorder is None:
            raise ValueError("Must pass -o with -f")
        from .handle_hd5f import main as handle_hdf5_main
        handle_hdf5_main(
            num_cores=args.num_cores,
            max_order=args.maxorder,
            max_homology_dim=args.maxdim,
            verbose=args.verbose,
            hdf5_folder=args.folder,
            chunksize=args.chunksize
        )

main()
