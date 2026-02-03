from .find_generating_subset import find_generating_subset
from mutable_lattice import Vector

def cover_submodule_with_actions(
    N: int,
    Zbasis: list[Vector],
    actions: list[Vector],
    e_to_Se: dict[int,tuple[int,...]],
    *,
    extra_greedy=False,
    ensure_minimal=False,
    verbose=False,
):
    def find_fixer_idempotent(vec):
        for e in e_to_Se:
            if vec.shuffled_by_action(actions[e]) == vec:
                return e
        raise AssertionError("Not a monoid?")
    id_to_idempotent = {id(vec): find_fixer_idempotent(vec) for vec in Zbasis}
    costs = [len(e_to_Se[id_to_idempotent[id(vec)]]) for vec in Zbasis]
    generating_subset = find_generating_subset(
        N, Zbasis, actions, costs,
        extra_greedy=extra_greedy,
        ensure_minimal=ensure_minimal,
        verbose=verbose)
    result_ZS_module = [id_to_idempotent[id(vec)] for vec in generating_subset]
    result_matrix = [
        vec.shuffled_by_action(actions[se])
        for vec, e in zip(generating_subset, result_ZS_module)
        for se in e_to_Se[e]
    ]
    return result_matrix, result_ZS_module

def make_actions(op, e_to_Se, e_to_s_to_ii, mod):
    N = sum(map(len, map(e_to_Se.get, mod)))
    n0 = 0
    _left_multiply_index_table = [[None]*N for s in range(len(op))]
    for e in mod:
        Se = e_to_Se[e]
        x_to_ii = e_to_s_to_ii[e]
        start_index_x_pairs = list(enumerate(Se, n0))
        for op_s, left_table_s in zip(op, _left_multiply_index_table):
            for start_index, x in start_index_x_pairs:
                left_table_s[start_index] = n0 + x_to_ii[op_s[x]]
        n0 += len(Se)
    return list(map(Vector, _left_multiply_index_table))

def cover_submodule(op, e_to_Se, e_to_x_to_ii, module, submodule_Zbasis, *, extra_greedy=False, ensure_minimal=False, verbose=False):
    N = sum(map(len, map(e_to_Se.get, module)))
    actions = make_actions(op, e_to_Se, e_to_x_to_ii, module)
    return cover_submodule_with_actions(
        N, submodule_Zbasis, actions, e_to_Se,
        extra_greedy=extra_greedy,
        ensure_minimal=ensure_minimal,
        verbose=verbose)
