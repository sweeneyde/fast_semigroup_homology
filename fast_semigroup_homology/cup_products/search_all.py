"""
Search for monoids with nontrivial cup products.
"""

from . import BarComplex
from .. import fast_integral_semigroup_homology

MAX_DIM = 4

def cup_product_exponents_worker(id_op):
    i, op = id_op
    # results.write(f"working on index {i} --> " + ";".join(''.join(map(str, row)) for row in op))
    homology_list = fast_integral_semigroup_homology(op, MAX_DIM)
    free_parts = [{0: h[0]} if 0 in h else {} for h in homology_list]
    torsion_parts = [{k: v for k, v in h.items() if k} for h in homology_list]
    cohomology_list = [t | f for t, f in zip([{}] + torsion_parts, free_parts)]
    bar = BarComplex(op)
    result = {}
    needed_pairs = [(a, b)
                    for a in range(1, MAX_DIM + 1)
                    for b in range(1, MAX_DIM + 1)
                    if a <= b and a + b <= MAX_DIM]
    needed_pairs.sort(key=sum)
    for a, b in needed_pairs:
        if not all(map(cohomology_list.__getitem__, [a, b, a + b])):
            continue
        a_invariants, b_invariants, exponents = bar.cup_product_exponents(a, b)
        if all(d == 1 for row in exponents for d in row):
            continue
        result[a, b] = exponents
    return i, cohomology_list, result

if __name__ == "__main__":
    from multiprocessing import Pool
    from pathlib import Path

    import h5py
    from tqdm import tqdm

    REPO = Path(__file__).parent.parent.parent
    SEMIGROUP_TABLES_FOLDER = (
        REPO.parent
        / "semisearch"
        / "results"
        / "monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units"
        # / "monoids_no_monoid_1sided_ideals_by_min_ideal_and_diagonal_and_units_bounded_qdiag"
    )

    with Pool(11) as pool:
        for order in [11]:
            nontrivial_results = {}
            with h5py.File(SEMIGROUP_TABLES_FOLDER / f"order{order}.hdf5") as tables_file:
                # Skip over groups: they're well-studied already.
                kind1_name, kind1_time, kind1_start, kind1_end = tables_file["kinds"][0]
                assert kind1_name == b"GROUPS"
                assert kind1_start == 0
                tables = tables_file["tables"]
                N = tables.shape[2]
                worker_data = ((i, tables[:,:,i]) for i in reversed(range(kind1_end, N)))
                results = pool.imap_unordered(cup_product_exponents_worker, worker_data)
                results = tqdm(results, f"order{order}", total=N-kind1_end, smoothing=0.0, miniters=1, mininterval=0.1)
                for i, cohomology_list, cups in results:
                    if cups:
                        nontrivial_results[i] = cups
                        results.write(f"NONTRIVIAL MATRIX: index {i}, {str(cohomology_list).replace(' ', '')}, {cups=}")
                        results.write("from " + ";".join(''.join(map(str, row)) for row in tables[:,:,i]))
            if nontrivial_results:
                print(f"{sorted(nontrivial_results)=}")
            else:
                print("no nontrivial results found")

# NONTRIVIAL MATRIX: index 6, [{0:1},{},{0:1},{},{2:1}], cups={(2, 2): [[2]]}
# from 010101;010110;232323;232332;012345;230154
# order6: 100%|████████████████████████████████████████████████████████████████| 5/5 [00:00<00:00, 18.69it/s]
# sorted(nontrivial_results)=[6]
# order7: 100%|██████████████████████████████████████████████████████████████| 24/24 [00:00<00:00, 91.34it/s]
# no nontrivial results found
# NONTRIVIAL MATRIX: index 243, [{0:1},{},{0:1},{},{4:1}], cups={(2, 2): [[2]]}
# from 01010101;01011010;23232323;23233232;01234567;23015674;01236745;23017456
# NONTRIVIAL MATRIX: index 246, [{0:1},{},{0:1},{},{2:2}], cups={(2, 2): [[2]]}
# from 01010011;01011100;23232233;23233322;01234567;01235476;23016745;23017654
# order8: 100%|███████████████████████████████████████████████████████████| 251/251 [00:01<00:00, 180.88it/s]
# sorted(nontrivial_results)=[243, 246]
# NONTRIVIAL MATRIX: index 3664, [{0:1},{},{0:1},{},{2:1}], cups={(2, 2): [[2]]}
# from 012301230;103210321;012310322;103201233;456745674;547654765;456754766;547645677;012345678
# NONTRIVIAL MATRIX: index 3663, [{0:1},{},{2:1,0:1},{},{2:1}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 012301230;103210321;012301232;103210323;456745674;547654765;456745676;547654767;012345678
# order9: 100%|██████████████████████████████████████████████████████████| 3663/3663 [00:42<00:00, 87.10it/s]
# sorted(nontrivial_results)=[3663, 3664]
# NONTRIVIAL MATRIX: index 66511, [{0:1},{},{0:1},{0:1},{2:1,0:2}], cups={(2, 2): [[2]]}
# from 0101010101;0101010110;2323232323;2323232332;0101010145;0101010154;2323232367;2323232376;0123456789;2301674598
# NONTRIVIAL MATRIX: index 67415, [{0:1},{},{0:1},{},{2:1}], cups={(2, 2): [[2]]}
# from 0101010101;0101010110;2323232323;2323232332;0101074147;2323632556;2323256365;0101410774;0123456789;2301674598
# NONTRIVIAL MATRIX: index 69542, [{0:1},{},{0:1},{},{2:1}], cups={(2, 2): [[2]]}
# from 0101001101;0101001110;2323223323;2323223332;0101446646;2323557757;0101446664;2323557775;0123456789;2301547698
# NONTRIVIAL MATRIX: index 69821, [{0:1},{},{0:1},{},{2:1}], cups={(2, 2): [[2]]}
# from 0101010101;0101101010;2323232323;2323323232;0123456789;2301547698;0123698547;2301789456;0123874965;2301965874
# NONTRIVIAL MATRIX: index 69824, [{0:1},{},{0:1},{},{6:1}], cups={(2, 2): [[2]]}
# from 0101010101;0101101010;2323232323;2323323232;0123456789;2301547698;0123678945;2301769854;0123894567;2301985476
# NONTRIVIAL MATRIX: index 71902, [{0:1},{},{2:1,0:1},{},{2:1}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012300;1032103211;0123012302;1032103213;4567456744;5476547655;4567456746;5476547657;0123012388;0123456789
# NONTRIVIAL MATRIX: index 71900, [{0:1},{},{2:1,0:1},{0:1},{2:1,0:1}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012300;1032103211;0123012302;1032103213;4567456744;5476547655;4567456746;5476547657;0123012308;0123456789
# NONTRIVIAL MATRIX: index 71901, [{0:1},{},{2:1,0:1},{0:1},{2:1,0:1}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012310;1032103201;0123012312;1032103203;4567456754;5476547645;4567456756;5476547647;1032103208;0123456789
# NONTRIVIAL MATRIX: index 71903, [{0:1},{},{2:1,0:1},{},{2:2}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012300;1032103211;0123012322;1032103233;4567456744;5476547655;4567456766;5476547677;0123456789;0123456798
# NONTRIVIAL MATRIX: index 71899, [{0:1},{},{0:2},{},{2:1}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012301;0123012310;0123012323;0123012332;4567456745;4567456754;4567456767;4567456776;0123456789;4567012398
# NONTRIVIAL MATRIX: index 71904, [{0:1},{},{2:1},{2:1},{2:1}], cups={(2, 2): [[2]]}
# from 0123012300;1032103211;0123012322;1032103233;4567456744;5476547655;4567456766;5476547677;0123456789;4567012398
# NONTRIVIAL MATRIX: index 71906, [{0:1},{},{2:1},{2:1},{2:1}], cups={(2, 2): [[2]]}
# from 0123012301;1032103210;0123012323;1032103232;4567456745;5476547654;4567456767;5476547676;0123456789;5476103298
# NONTRIVIAL MATRIX: index 71905, [{0:1},{},{2:1,0:1},{},{2:2}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012301;1032103210;0123012323;1032103232;4567456745;5476547654;4567456767;5476547676;0123456789;1032547698
# NONTRIVIAL MATRIX: index 71907, [{0:1},{},{2:1,0:1},{},{2:2}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012302;1032103213;0123012320;1032103231;4567456746;5476547657;4567456764;5476547675;0123456789;4567012398
# NONTRIVIAL MATRIX: index 71910, [{0:1},{},{0:1},{0:1},{2:1,0:1}], cups={(2, 2): [[2]]}
# from 0123012310;1032103201;0123103212;1032012303;4567456754;5476547645;4567547656;5476456747;1032103208;0123456789
# NONTRIVIAL MATRIX: index 71911, [{0:1},{},{0:1},{},{2:1}], cups={(2, 2): [[2]]}
# from 0123012300;1032103211;0123103202;1032012313;4567456744;5476547655;4567547646;5476456757;0123012388;0123456789
# NONTRIVIAL MATRIX: index 71908, [{0:1},{},{2:1,0:1},{},{2:2}], cups={(2, 2): [[2, 2], [2, 2]]}
# from 0123012303;1032103212;0123012321;1032103230;4567456747;5476547656;4567456765;5476547674;0123456789;5476103298
# NONTRIVIAL MATRIX: index 71912, [{0:1},{},{0:1},{},{2:2}], cups={(2, 2): [[2]]}
# from 0123012300;1032103211;0123103222;1032012333;4567456744;5476547655;4567547666;5476456777;0123456789;0123456798
# NONTRIVIAL MATRIX: index 71909, [{0:1},{},{0:1},{0:1},{2:1,0:1}], cups={(2, 2): [[2]]}
# from 0123012300;1032103211;0123103202;1032012313;4567456744;5476547655;4567547646;5476456757;0123012308;0123456789
# NONTRIVIAL MATRIX: index 71914, [{0:1},{},{0:1},{},{4:1}], cups={(2, 2): [[4]]}
# from 0123012301;1032103210;0123103223;1032012332;4567456745;5476547654;4567547667;5476456776;0123456789;1032547698
# order10: 100%|███████████████████████████████████████████████████████| 71914/71914 [27:17<00:00, 43.91it/s]
# sorted(nontrivial_results)=[66511, 67415, 69542, 69821, 69824, 71899, 71900, 71901, 71902, 71903, 71904, 71905, 71906, 71907, 71908, 71909, 71910, 71911, 71912, 71914]