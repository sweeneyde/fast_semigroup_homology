from pathlib import Path
import re
from contextlib import nullcontext
from collections import Counter
from functools import cache

import ast
from tqdm import tqdm
import numpy as np

from .homology import fast_integral_semigroup_homology

@cache
def str_trivial_homology(maxdim):
    return str([{0: 1}] + [{}] * maxdim)

def monoid_homology_worker(index_table_maxdim_verbose):
    index, op, maxdim, verbose = index_table_maxdim_verbose
    homology_list = fast_integral_semigroup_homology(
        op, maxdim=maxdim, verbose=verbose)
    return index, str(homology_list)

def hdf5_compute_homology(*, multiprocessing_mapper, input_folder, max_order, max_homology_dim, verbose):
    import h5py # import here so the multiprocessing workers never need to import h5py
    input_folder = Path(input_folder)
    integer_to_path = {}
    for filepath in input_folder.iterdir():
        order = int(re.fullmatch(r"order(\d+)\.hdf5", filepath.name).group(1))
        if order <= max_order:
            integer_to_path[order] = filepath
    for order in range(1, max_order + 1):
        assert order in integer_to_path, (order, integer_to_path)
    integer_to_path = dict(sorted(integer_to_path.items()))
    output_folder = Path(__file__).parent.parent / "results" / input_folder.name
    output_folder.mkdir(exist_ok=True)
    output_filename = f"maxorder{max_order}_maxdim{max_homology_dim}.hdf5"
    output_path = output_folder / output_filename
    with h5py.File(output_path, "x") as output_hfile:
        hgl_dset = output_hfile.create_dataset(
            name="homology_group_lists",
            shape=(1,),
            maxshape=(None,),
            dtype=h5py.string_dtype(),
        )
        hgl_dset[0] = "(not yet computed)"
        hgl_to_index = {}
        def get_hgl_index(homology_group_list):
            index = hgl_to_index.get(homology_group_list)
            if index is not None:
                return index
            result = len(hgl_dset)
            hgl_dset.resize((result+1,))
            hgl_dset[result] = homology_group_list
            hgl_to_index[homology_group_list] = result
            return result

        for order, filepath in integer_to_path.items():
            print(f"working on {filepath.name}...")
            with h5py.File(filepath, "r") as input_hfile:
                input_tables = input_hfile["tables"]
                input_kinds = input_hfile["kinds"]
                shape = input_tables.shape
                assert len(shape) == 3
                assert shape[0] == shape[1] == order
                count = shape[2]

                output_dset = output_hfile.create_dataset(
                    name=f"order{order}",
                    shape=(count,),
                    dtype="u4",
                    compression="gzip",
                    compression_opts=9,
                    fletcher32=True,
                    shuffle=True,
                )
                progress_bar = tqdm(desc=filepath.name, total=count, miniters=1, smoothing=0.0)
                def tables_iterator():
                    for name, _, i0, i1 in input_kinds:
                        name = name.decode("ascii")
                        if name.startswith("min_1_"):
                            assert i0 < i1
                            table = input_tables[:,:,i0].tolist()
                            _, string = monoid_homology_worker((i0, table, max_homology_dim, verbose))
                            output_dset[i0:i1] = get_hgl_index(string)
                            progress_bar.update(i1-i0)
                        else:
                            for i in range(i0, i1):
                                yield i, input_tables[:,:,i].tolist(), max_homology_dim, verbose
                worker = monoid_homology_worker
                finished_work = multiprocessing_mapper(worker, tables_iterator())
                for index, string in finished_work:
                    if verbose:
                        print(f"===========================")
                        print(f"==== {order}_{index}: {string}")
                        print(f"===========================")
                    output_dset[index] = get_hgl_index(string)
                    progress_bar.update(1)
                progress_bar.close()
            print(f"done with {filepath.name}")
    return output_path

def summarize_hdf5_as_markdown(hfile_path, max_homology_dim, dt, num_cores):
    SUPER_DIGITS = "⁰¹²³⁴⁵⁶⁷⁸⁹"
    SUB_DIGITS = "₀₁₂₃₄₅₆₇₈₉"
    @cache
    def superscript(n):
        return "".join(SUPER_DIGITS[int(c)] for c in str(n))
    def maybe_superscript(n):
        assert n >= 1
        return "" if n == 1 else superscript(n)
    @cache
    def subscript(n):
        return "".join(SUB_DIGITS[int(c)] for c in str(n))
    def H_sub(n):
        return "𝐻" + subscript(n)
    @cache
    def summand_name(n, count):
        if n == 0:
            return "ℤ" + maybe_superscript(count)
        else:
            return "𝐶" + subscript(n) + maybe_superscript(count)
    def group_name(group):
        if not group:
            return "·"
        return "×".join(summand_name(n, count) for n, count in group)

    import h5py
    with h5py.File(hfile_path, "r") as hfile:
        dset_names = set(hfile.keys())
        dset_names.remove("homology_group_lists")
        homology_group_lists_dset = hfile["homology_group_lists"]
        dset_names = sorted(dset_names, key=lambda name: int(name.removeprefix("order")))
        md_filepath = hfile_path.parent / (hfile_path.name.removesuffix(".hdf5") + ".md")
        with open(md_filepath, "w", encoding="utf-8") as f:
            print("# Lists of Homology Groups", file=f)
            print(f"Computation wall time with {num_cores} cores: `{dt}`", file=f)
            for dset_name in dset_names:
                dset = hfile[dset_name]
                index_counts = Counter()
                progress_bar = tqdm(desc="counting", smoothing=0.0, miniters=1, total=len(dset))
                for _slice in dset.iter_chunks():
                    arr = dset[_slice]
                    unique, counts = np.unique(arr, return_counts=True)
                    for u, c in zip(unique, counts):
                        index_counts[u] += c
                    progress_bar.update(len(arr))
                progress_bar.close()
                print(f"\n## {dset_name}", file=f)
                print(file=f)
                print(*(["Count"] + [H_sub(i) for i in range(max_homology_dim + 1)]), sep=" | ", file=f)
                print(*(["--:"] + [":--:" for i in range(max_homology_dim + 1)]), sep=" | ", file=f)
                list_counts = []
                for index, count in index_counts.items():
                    assert count > 0
                    homology_group_list_str = homology_group_lists_dset[index].decode("ascii")
                    homology_group_list = ast.literal_eval(homology_group_list_str)
                    list_counts.append(([list(d.items()) for d in homology_group_list], count))
                for homology_group_list, count in sorted(list_counts):
                    row = [f"{count:,}"] + list(map(group_name, homology_group_list))
                    print(*row, sep=" | ", file=f)
    print(f"Wrote to {md_filepath}.")

def main(*, num_cores, max_order, max_homology_dim, verbose, hdf5_folder, chunksize):
    from time import time
    t0 = time()
    if num_cores > 1:
        import multiprocessing as mp
        mp.set_start_method("spawn")
        context = mp.Pool(num_cores)
    else:
        context = nullcontext()
    with context as pool:
        if num_cores > 1:
            def mapper(func, it):
                return pool.imap_unordered(func, it, chunksize=chunksize)
        else:
            mapper = map
        output_path = hdf5_compute_homology(
            multiprocessing_mapper=mapper,
            input_folder=hdf5_folder,
            max_order=max_order,
            max_homology_dim=max_homology_dim,
            verbose=verbose,
        )
    t1 = time()
    from datetime import timedelta
    dt = str(timedelta(seconds=t1-t0))
    print(f"{dt} elapsed")
    summarize_hdf5_as_markdown(output_path, max_homology_dim, dt, num_cores)
