"""
Search for monoids with nontrivial rational cohomology cup products.

So far, I have found no rational cup products:
    * H1 = abelianization(groupCompletion(M)) is finite; no free parts
    * from monoids of order <= 11, no H2 x H2 --> H4.
    * from monoids of order <= 10, no H2 x H3 --> H5.
    * from monoids of order <= 9, no H3 x H3 --> H6 nor H2 x H4 --> H6 nor H3 x H4 --> H7
    * from monoids of order <= 7, no H2 x H5 --> H7 nor H3 x H5 --> H8 nor H4 x H5 --> H9 nor H5 x H5 --> H10
    * from monoids of order <= 15 with bounded_qdiag, from H2 x H2 --> H4

The first nontrivial cup product I know of is from the 25-element monoid
    Rect(2,2)^1 x Rect(2,2)^1, with classifying space ~= S2 x S2.
"""

from . import cup_product_matrix

DIM_A = 2
DIM_B = 2

def cup_product_matrix_worker(id_op):
    i, op = id_op
    return i, cup_product_matrix(op, DIM_A, DIM_B)

if __name__ == "__main__":
    from ast import literal_eval
    from multiprocessing import Pool
    from pathlib import Path
    import random

    import h5py
    from tqdm import tqdm

    REPO = Path(__file__).parent.parent.parent
    HOMOLOGY_RESULTS = (
        REPO
        / "results"
        # / "monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units_bounded_qdiag"
        # / "maxorder15_maxdim4.hdf5"
        / "monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units"
        / "refined_maxorder11_maxdim4.hdf5"
    )
    SEMIGROUP_TABLES_FOLDER = (
        REPO.parent
        / "semisearch"
        / "results"
        / "monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units"
        # / "monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units_bounded_qdiag"
    )

    MAXORDER = int(HOMOLOGY_RESULTS.stem.rpartition("_")[0].partition("maxorder")[2])

    with Pool(11) as pool:
        with h5py.File(HOMOLOGY_RESULTS) as hr_file:
            good_codes = {
                code for code, s in enumerate(hr_file["homology_group_lists"])
                if code > 0
                and (h := literal_eval(s.decode("ascii")))
                and 0 in h[DIM_A]
                and 0 in h[DIM_B]
                and 0 in h[DIM_A + DIM_B]
            }
            for order in range(1, MAXORDER + 1):
                nontrivial_results = {}
                hr_dset = hr_file[f"order{order}"]
                good_ids = [i for i, code in enumerate(hr_dset) if code in good_codes]
                random.Random(0).shuffle(good_ids)
                with h5py.File(SEMIGROUP_TABLES_FOLDER / f"order{order}.hdf5") as tables_file:
                    tables = tables_file["tables"]
                    if tables.shape[2] != len(hr_dset):
                        raise ValueError("inconsistent input data files")
                    worker_data = ((i, tables[:,:,i]) for i in good_ids)
                    results = pool.imap_unordered(cup_product_matrix_worker, worker_data)
                    results = tqdm(results, f"order{order}", total=len(good_ids), smoothing=0.0, miniters=1, mininterval=0.1)
                    for i, matrix in results:
                        if any(int(z) for x in matrix for y in x for z in y):
                            nontrivial_results[i] = matrix
                            print(f"FOUND NONTRIVIAL MATRIX: {matrix} for index {i}")
                if nontrivial_results:
                    print(f"{nontrivial_results=}")
                else:
                    print("no nontrivial results found")
