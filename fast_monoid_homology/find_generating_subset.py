from mutable_lattice import Vector, Lattice

def find_generating_subset(
    Zbasis: list[Vector],
    actions: list[Vector],
    costs: list[int],
    extra_greedy=False,
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

    assert set(map(len, Zbasis)) == {N}
    assert set(map(len, actions)) == {N}
    for act in actions:
        assert set(act) <= set(range(N))
    assert len(costs) == R

    if not extra_greedy:
        def sort_by_increasing_cost():
            id_to_cost = {id(vec): c for vec, c in zip(Zbasis, costs)}
            return sorted(Zbasis, key=lambda vec: id_to_cost[id(vec)])
        Zbasis = sort_by_increasing_cost()
        if verbose:
            print(f"sorted by increasing cost")

    def relativize():
        """
        Convert to the provided basis,
        so the vectors we use are smaller,
        and so there will be fewer columns without a pivot.
        """
        K = Lattice(N, Zbasis)
        assert K.rank == R
        relative_lattices = []
        relative_vectors = []
        for vec in Zbasis:
            relative_vectors.append(K.coefficients_of(vec))
            translated = [K.coefficients_of(vec.shuffled_by_action(act)) for act in actions]
            rel_lattice = Lattice(R, translated, maxrank=len(actions))
            # trim to save a tiny amount of memory
            rel_lattice = Lattice(R, rel_lattice.get_basis(), maxrank=rel_lattice.rank)
            relative_lattices.append(rel_lattice)
        return relative_vectors, relative_lattices
    relative_vectors, relative_lattices = relativize()
    if verbose:
        print(f"relativized.")

    def shuffle_columns():
        """
        Shuffle all vectors in the same way,
        so that the least used columns come first.
        """
        nonlocal relative_vectors, relative_lattices
        column_to_count = [0] * R
        for rel_lattice in relative_lattices:
            for rel_vec in rel_lattice.get_basis():
                for i in filter(rel_vec.__getitem__, range(R)):
                    column_to_count[i] += 10_000_000 + abs(rel_vec[i]).bit_length()
        columns_by_increasing_count = sorted(range(R), key=column_to_count.__getitem__)
        sort_action = [None] * R
        for i, col in enumerate(columns_by_increasing_count):
            sort_action[col] = i
        sort_action = Vector(sort_action)
        relative_vectors = [v.shuffled_by_action(sort_action) for v in relative_vectors]
        relative_lattices = [Lattice(R, [v.shuffled_by_action(sort_action)
                                        for v in reversed(rel_lattice.get_basis())],
                                    maxrank=rel_lattice.rank)
                             for rel_lattice in relative_lattices]
    shuffle_columns()
    if verbose:
        print(f"shuffled columns.")

    def one_pass_cover(existing_solution):
        """
        Walk over the list of vectors, keeping only the vectors
        not already covered by previous vectors.
        """
        new_solution = []
        L = Lattice(R)
        for i in existing_solution:
            if relative_vectors[i] not in L:
                L += relative_lattices[i]
                new_solution.append(i)
            else:
                relative_lattices[i] = None
        if verbose:
            print(f"shrank solution {len(existing_solution)} --> {len(new_solution)}")
        return new_solution

    def greedy_cover():
        """Add in the most efficient vector at every step"""
        solution = []
        uncovered = set(range(R))
        L = Lattice(R)
        while uncovered:
            best_efficiency = -1
            best_newly_covered = None
            best_index = None
            best_hypothetical_sum = None
            for i in uncovered:
                hypothetical_sum = L + relative_lattices[i]
                newly_covered = [j for j in uncovered if relative_vectors[j] in hypothetical_sum]
                efficiency = len(newly_covered) / (0.001 + costs[i])
                if efficiency > best_efficiency:
                    best_efficiency = efficiency
                    best_newly_covered = newly_covered
                    best_hypothetical_sum = hypothetical_sum
                    best_index = i
            solution.append(best_index)
            uncovered.difference_update(best_newly_covered)
            L = best_hypothetical_sum
        if verbose:
            print(f"Greedy solution has {len(solution)} vectors")
        return solution

    solution = greedy_cover() if extra_greedy else list(range(R))
    solution = one_pass_cover(solution)
    solution = one_pass_cover(solution[::-1])
    # (could shuffle the solution and do a few more one-pass covers here)

    def do_ensure_minimal(existing_solution):
        """
        Keep only those vectors not covered by the rest
        of the vectors in the solution.
        """
        new_solution = []
        suffix_sums = [Lattice(R, maxrank=0)]
        for i in reversed(existing_solution):
            suffix_sums.append(suffix_sums[-1] + relative_lattices[i])
        prefix_sum = Lattice(R)
        for j, i in enumerate(existing_solution):
            rel_vec = relative_vectors[i]
            suff = suffix_sums[len(existing_solution) - j - 1]
            if rel_vec in (prefix_sum + suff):
                # redundant
                pass
            else:
                prefix_sum += relative_lattices[i]
                new_solution.append(i)
        if verbose:
            print(f"minimality pass: {len(existing_solution)} --> {len(new_solution)}")
        return new_solution
    if ensure_minimal:
        solution = do_ensure_minimal(solution)
    solution.sort()
    return [Zbasis[i] for i in solution]