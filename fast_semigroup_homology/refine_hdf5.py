import ast
from collections import Counter, defaultdict
import shutil

import numpy as np
import h5py
from mutable_lattice import relations_among

from .handle_hd5f import summarize_hdf5_as_markdown
from .homology import integral_monoid_homology

class KernelJobTooBig(Exception):
    pass

def kernel_bounded(vectors, bounds):
    if not vectors:
        return []
    R = len(vectors)
    N = len(vectors[0])
    bits = max(max(map(abs, v)) for v in vectors).bit_length()
    max_N, max_R, max_bits = bounds
    if R > max_R or N > max_N or bits > max_bits:
        raise KernelJobTooBig(f"{R=}{'(over)' if R > max_R else ''},"
                              f"{N=}{'(over)' if N > max_N else ''},"
                              f"{bits=}{'(over)' if bits > max_bits else ''}")
    return relations_among(vectors).get_basis()

def homology_worker_small_kernel(index_op_maxdim_bounds):
    index, op, maxdim, bounds = index_op_maxdim_bounds
    try:
        def my_kernel_implementation(vectors, verbose=False):
            return kernel_bounded(vectors, bounds)
        hgl = integral_monoid_homology(
            op, maxdim,
            kernel_implementation=my_kernel_implementation,
        )
        return index, str(hgl)
    except KernelJobTooBig as e:
        print(f"index {index}:", e)
        return index, "KernelJobTooBig"

def attempt_job_in_parallel(*,
        mapper,
        work_iterator,
):
    extended_hgl_to_semigroup_index_list = defaultdict(list)
    for index, result in mapper(homology_worker_small_kernel, work_iterator):
        if result == "KernelJobTooBig":
            return None
        else:
            extended_hgl_to_semigroup_index_list[result].append(index)
    for semigroup_index_list in extended_hgl_to_semigroup_index_list.values():
        semigroup_index_list.sort()
    return extended_hgl_to_semigroup_index_list

def hgl_shorthand(homology_group_list):
    # Just take up less space in the print output
    return ",".join("+".join(f"{k if k else 'Z'}" + (f"^{v}" if v > 1 else "")
                             for k, v in gp.items())
                    for gp in homology_group_list
    ).join("[]")

def hdf5_refine_homology(*,
        num_cores,
        input_folder,
        existing_homology_filepath,
        maxdim_limit,
        max_R,
        max_N,
        max_bits,
        max_count_to_refine,
):
    import multiprocessing as mp
    mp.set_start_method("spawn")
    pool = mp.Pool(num_cores)
    bounds = max_N, max_R, max_bits
    output_filepath = existing_homology_filepath.parent / ("refined_" + existing_homology_filepath.name)
    if output_filepath.exists():
        raise ValueError(f"already exists: {output_filepath}")
    shutil.copy(existing_homology_filepath, output_filepath)
    del existing_homology_filepath
    with h5py.File(output_filepath, "r+") as output_hfile:
        hgl_dset = output_hfile["homology_group_lists"]
        dset_names = set(output_hfile.keys())
        orders = sorted(int(name.removeprefix("order"))
                            for name in dset_names
                            if name != "homology_group_lists")
        hgl_str_list = [hgl.decode("ascii") for hgl in hgl_dset]
        hgl_lists = [None] + list(map(ast.literal_eval, hgl_str_list[1:]))
        hgl_str_to_index = {hgl: i for i, hgl in enumerate(hgl_str_list)}
        def get_hgl_index(homology_group_list):
            assert isinstance(homology_group_list, str)
            index = hgl_str_to_index.get(homology_group_list)
            if index is not None:
                return index
            result = len(hgl_dset)
            hgl_dset.resize((result+1,))
            hgl_str_list.append(homology_group_list)
            hgl_lists.append(ast.literal_eval(homology_group_list))
            hgl_dset[result] = homology_group_list
            hgl_str_to_index[homology_group_list] = result
            return result
        for order in orders:
            print(f"\n\n===== order{order} =====")
            print("First scan to count occurrences of homology lists...")
            hgl_index_dset = output_hfile[f"order{order}"]
            index_counts = Counter()
            for (_slice,) in hgl_index_dset.iter_chunks():
                for u, c in zip(*np.unique(hgl_index_dset[_slice], return_counts=True)):
                    index_counts[u] += c

            # Stack of work to do: each stack item is
            #   count, hgl_index, (either None or a semigroup index list)
            stack = [[count, int(hgl_index), None]
                     for hgl_index, count in index_counts.items()
                     if count <= max_count_to_refine]
            stack.sort(key=lambda stackitem: stackitem[0], reverse=True)

            print("Second scan find indices of rare homology lists...")
            rare_hgl_index_to_semigroup_indexes = {}
            for stackitem in stack:
                count, hgl_index, _ = stackitem
                if count <= 1000:
                    stackitem[2] = rare_hgl_index_to_semigroup_indexes[hgl_index] = []
            rare_hgl_indexes = sorted(rare_hgl_index_to_semigroup_indexes.keys())
            for (_slice,) in hgl_index_dset.iter_chunks():
                chunk = hgl_index_dset[_slice]
                (rare,) = np.isin(chunk, rare_hgl_indexes).nonzero()
                for i in rare:
                    rare_hgl_index_to_semigroup_indexes[chunk[i]].append(int(_slice.start + i))
            del rare_hgl_index_to_semigroup_indexes # the lists are still on the stack.

            tables_file_path = input_folder / f"order{order}.hdf5"
            with h5py.File(tables_file_path, "r") as tables_file:
                tables = tables_file["tables"]
                if tables.shape[2] != len(hgl_index_dset):
                    raise ValueError(
                        f"{tables_file_path} had {tables.shape[2]} tables (shape={tables.shape}), but "
                        f"{f"order{order}"} dset had {len(hgl_index_dset)}"
                    )
                while stack:
                    print(f"{len(stack)=}")
                    count, hgl_index, semigroup_indexes = stack.pop()
                    existing_hgl = hgl_lists[hgl_index]
                    if len(existing_hgl) - 1 >= maxdim_limit:
                        print(f"{existing_hgl} already has maxdim>={maxdim_limit}")
                        continue
                    print(f"Attempting to extend {hgl_shorthand(existing_hgl)} x{count}...")
                    if semigroup_indexes is None:
                        semigroup_indexes = []
                        for (_slice,) in hgl_index_dset.iter_chunks():
                            chunk = hgl_index_dset[_slice]
                            (matching,) = (chunk == hgl_index).nonzero()
                            semigroup_indexes.extend(_slice.start + matching)
                    assert len(semigroup_indexes) == count
                    # Go one further--this length includes dimension 0
                    maxdim = len(existing_hgl)
                    def work_iterator():
                        for i in semigroup_indexes:
                            yield i, tables[:,:,i], maxdim, bounds
                    extended_hgl_to_semigroup_index_list = attempt_job_in_parallel(
                        mapper=pool.imap_unordered if len(semigroup_indexes) > 1 else map,
                        work_iterator=work_iterator(),
                    )
                    if extended_hgl_to_semigroup_index_list is None:
                        print(f"abandoned work on {hgl_shorthand(existing_hgl)}")
                        # abandon this pool and make a new one
                        if len(semigroup_indexes) > 1:
                            pool.terminate()
                            pool = mp.Pool(num_cores)
                        continue

                    all_returned_semigroup_indexes = []
                    for sublist in extended_hgl_to_semigroup_index_list.values():
                        all_returned_semigroup_indexes.extend(sublist)
                    all_returned_semigroup_indexes.sort()
                    if all_returned_semigroup_indexes != semigroup_indexes:
                        raise ValueError(f"Job did not partition semigroup_indexes")
                    del all_returned_semigroup_indexes

                    results = sorted(extended_hgl_to_semigroup_index_list.items(),
                                    key=lambda x: len(x[1]),
                                    reverse=True)
                    if len(results) >= 2:
                        print("NEW SPLIT:")
                    for extended_hgl_str, semigroup_indexes in results:
                        extended_hgl = ast.literal_eval(extended_hgl_str)
                        assert extended_hgl[:len(existing_hgl)] == existing_hgl
                        print(f"found {hgl_shorthand(extended_hgl)} x{len(semigroup_indexes)}")
                        hgli = get_hgl_index(str(extended_hgl))
                        hgl_index_dset[semigroup_indexes] = hgli
                        if len(extended_hgl)-1 < maxdim_limit:
                            # Extend the stack to try to go further next time.
                            stack.append([len(semigroup_indexes), hgli, semigroup_indexes])
    pool.terminate()
    return output_filepath

def trim_unused_homology_group_lists(filepath):
    print("removing unused indexes...")
    copy_filepath = filepath.parent / ("copy_" + filepath.name)
    with (
        h5py.File(copy_filepath, "x") as hfile_copy,
        h5py.File(filepath, "r") as hfile,
    ):
        print(f"finding used indexes...")
        used = {0} # always include the "not yet computed" bin
        for name, hgl_index_dset in hfile.items():
            if name != "homology_group_lists":
                for (_slice,) in hgl_index_dset.iter_chunks():
                    used.update(np.unique(hgl_index_dset[_slice]))
        hgl_dset = hfile["homology_group_lists"]
        print(f"shrinking {len(hgl_dset)} --> {len(used)}")
        hgl_dset_copy = hfile_copy.create_dataset_like(
            "homology_group_lists",
            hgl_dset,
            shape=(len(used),),
        )
        print("copying homology_group_lists")
        new_to_old = sorted(used)
        for new, old in enumerate(new_to_old):
            hgl_dset_copy[new] = hgl_dset[old]
        old_to_new = np.zeros(shape=(len(hgl_dset),))
        for new, old in enumerate(new_to_old):
            old_to_new[old] = new
        for name, hgl_index_dset in hfile.items():
            if name == "homology_group_lists":
                continue
            print(f"copying {name}")
            hgl_index_dset_copy = hfile_copy.create_dataset_like(
                name,
                hgl_index_dset,
            )
            for (_slice,) in hgl_index_dset.iter_chunks():
                hgl_index_dset_copy[_slice] = old_to_new[hgl_index_dset[_slice]]
    print("Transferring back...")
    filepath.unlink()
    copy_filepath.rename(filepath)
    print("Done trimming.")

def main(*,
        input_folder,
        existing_homology_filepath,
        num_cores,
        max_homology_dim,
        max_R,
        max_N,
        max_bits,
        max_count_to_refine,
):
    from time import time
    t0 = time()
    output_path = hdf5_refine_homology(
        num_cores=num_cores,
        input_folder=input_folder,
        existing_homology_filepath=existing_homology_filepath,
        maxdim_limit=max_homology_dim,
        max_R=max_R,
        max_N=max_N,
        max_bits=max_bits,
        max_count_to_refine=max_count_to_refine,
    )
    trim_unused_homology_group_lists(output_path)
    t1 = time()
    from datetime import timedelta
    dt = str(timedelta(seconds=t1-t0))
    print(f"{dt} elapsed")
    info = {
        "max_homology_dim": max_homology_dim,
        "max_R": max_R,
        "max_N": max_N,
        "max_bits": max_bits,
        "max_count_to_refine": max_count_to_refine,
    }
    extra_info_lines = ["Refined using parameters:"] + [
        f"* {key}: {val}" for key, val in info.items()
    ]
    summarize_hdf5_as_markdown(
        output_path,
        max_homology_dim,
        dt,
        num_cores,
        extra_info_lines=extra_info_lines
    )
