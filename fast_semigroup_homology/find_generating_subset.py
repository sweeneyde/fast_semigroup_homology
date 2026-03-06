from mutable_lattice import Vector, Lattice
from functools import cache

def find_generating_subset(
    Zbasis: list[Vector],
    actions: list[Vector],
    costs: list[int],
    ensure_minimal=False,
    verbose=False,
):
    """
    Given a Z-basis for a sublattice of Z^N, along with
    a list of "shuffling action functions" to apply to the vectors,
    find a subset of the vectors such that the whole sublattice
    is spanned by the subset and its images under the actions.

    Zbasis is a list of Vectors of length N.
    actions is a list of Vectors of length N with entries in range(N),
        describing a collection of ways to shuffle around the
        entries of the Zbasis. The result of the shuffles should
        remain in the span of the Zbasis.
    costs is a list of integers: we aim to minimize
        the total cost of the returned subset.
    ensure_minimal is a bool: if True, we guarantee that our subset
        is minimal with respect to inclusion.
    """
    R = len(Zbasis)
    if R == 0:
        return []
    N = len(Zbasis[0])
    if verbose:
        print(f"Covering a rank-{R} sublattice of Z^{N}")

    # assert set(map(len, Zbasis)) == {N}
    # assert set(map(len, actions)) == {N}
    # for act in actions:
    #     assert set(act) <= set(range(N))
    # assert len(costs) == R

    def sort_by_increasing_cost():
        id_to_cost = {id(vec): c for vec, c in zip(Zbasis, costs)}
        return sorted(Zbasis, key=lambda vec: id_to_cost[id(vec)])
    Zbasis = sort_by_increasing_cost()
    if verbose:
        print(f"sorted by increasing cost.")

    K = Lattice(N, Zbasis)
    relative_vectors = [K.coefficients_of(vec) for vec in Zbasis]
    @cache
    def relative_lattice(i):
        vec = Zbasis[i]
        translated = [K.coefficients_of(vec.shuffled_by_action(act))
                      for act in actions]
        return Lattice(R, translated, maxrank=len(actions))

    def one_pass_cover(existing_solution):
        """
        Walk over the list of vectors, keeping only the vectors
        not already covered by previous vectors.
        """
        new_solution = []
        L = Lattice(R)
        it = existing_solution
        if verbose:
            from tqdm import tqdm
            it = tqdm(existing_solution, "cover", miniters=1)
        for i in it:
            if relative_vectors[i] not in L:
                L += relative_lattice(i)
                new_solution.append(i)
        if verbose:
            print(f"shrank solution {len(existing_solution)} --> {len(new_solution)}")
        return new_solution

    solution = one_pass_cover(list(range(R)))
    solution = one_pass_cover(solution[::-1])

    def do_ensure_minimal(existing_solution):
        """
        Keep only those vectors not covered by the rest
        of the vectors in the solution.
        """
        new_solution = []
        suffix_sums = [Lattice(R, maxrank=0)]
        for i in reversed(existing_solution):
            suffix_sums.append(suffix_sums[-1] + relative_lattice(i))
        prefix_sum = Lattice(R)
        it = existing_solution
        if verbose:
            from tqdm import tqdm
            it = tqdm(existing_solution, "minimality", miniters=1)
        for j, i in enumerate(it):
            rel_vec = relative_vectors[i]
            suff = suffix_sums[len(existing_solution) - j - 1]
            if rel_vec not in (prefix_sum + suff):
                new_solution.append(i)
                prefix_sum += relative_lattice(i)
        if verbose:
            print(f"minimality pass: {len(existing_solution)} --> {len(new_solution)}")
        return new_solution
    if ensure_minimal:
        solution = do_ensure_minimal(solution)
    solution.sort()
    return [Zbasis[i] for i in solution]