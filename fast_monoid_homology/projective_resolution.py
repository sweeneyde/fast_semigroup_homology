from collections import Counter
from functools import cache
from .cover_submodule import (
    cover_submodule_with_actions,
)
from .kernels import (
    mutable_lattice_kernel,
    mutable_lattice_kernel_with_reshuffling,
    mutable_lattice_kernel_with_pre_hnf,
    sage_kernel,
)
from .normalized_invariants import invariant_factors
from mutable_lattice import Vector, Lattice

class ProjectiveResolution:
    """Data representing a projective resolution of ZZ[S]-modules for some monoid S"""
    __slots__ = [
        # the monoid multiplication table as a tuple[tuple[int]]
        "op",
        # the integer index serving as an identity element for the monoid
        "identity",
        # The multiplication table (S elements) x (Y elements) --> (Y elements), as list[Vector]
        "left_S_set_action",
        # a mapping int -> tuple[int] from idempotents e to the tuple of distinct elements op[...][e]
        "e_to_Se",
        # int -> (int -> int). This stores the index at which element in Se appears in Se, for each e.
        "e_to_s_to_ii",
        # A mapping from ResolutionNode arguments to ResolutionNode, allowing re-use
        "node_cache",
        # The ResolutionNode for dimension 0 of this resolution
        "root",
        # A function that finds the relations among a given list of Vectors.
        "kernel_implementation",
    ]

    def __init__(self, op, *, left_S_set_action=None, check=True, kernel_implementation=None):
        if left_S_set_action is None:
            left_S_set_action = [[0] for _ in range(len(op))]
        left_S_set_action = [Vector(list(act)) for act in left_S_set_action]
        op = tuple(map(tuple, op))
        S = range(len(op))
        idempotents = [e for e in S if op[e][e] == e]
        [identity] = [e for e in idempotents if all(op[x][e] == x == op[e][x] for x in S)]
        Y = range(len(left_S_set_action[0]))
        if check:
            if len(left_S_set_action) != len(op):
                raise ValueError(f"{len(left_S_set_action)=} mismatches {len(op)=}")
            if set(map(len, op)) != {len(op)}:
                raise ValueError("op was not square")
            if len(set(map(len, left_S_set_action))) != 1:
                raise ValueError("inconsistent numbers of states")
            act = left_S_set_action
            for i in S:
                for j in S:
                    ij = op[i][j]
                    for k in S:
                        jk = op[j][k]
                        if op[ij][k] != op[ij][k]:
                            raise ValueError(f"op[op[{i}][{j}]][{k}]=op[{ij}][{k}]={op[ij][k]}, but "
                                             f"op[{i}][op[{j}][{k}]]=op[{i}][{jk}]={op[i][jk]}")
                    for y in Y:
                        jy = act[j][y]
                        if act[ij][y] != act[i][jy]:
                            raise ValueError(f"act[op[{i}][{j}]][{k}]=act[{ij}][{k}]={act[ij][k]}, but "
                                             f"act[{i}][act[{j}][{k}]]=act[{i}][{jy}]={act[i][jy]}")
            if act[identity] != Vector(list(Y)):
                raise ValueError(f"act[identity]=act[{identity}]={act[identity]} was not identity")
        e_Se_pairs = [(e, sorted({op[x][e] for x in S})) for e in idempotents]
        e_to_Se = dict(sorted(e_Se_pairs, key=lambda e_Se: len(e_Se[1])))
        e_to_s_to_ii = {
            e: {x: ii for ii, x in enumerate(Se)}
            for e, Se in e_to_Se.items()
        }
        self.op = op
        self.identity = identity
        self.left_S_set_action = left_S_set_action
        self.e_to_Se = e_to_Se
        self.e_to_s_to_ii = e_to_s_to_ii
        if kernel_implementation is None:
            kernel_implementation = mutable_lattice_kernel
        self.kernel_implementation = kernel_implementation

        augmentation_module_Zbasis = Lattice.full(len(Y)).get_basis()
        mat0, mod0 = cover_submodule_with_actions(
            augmentation_module_Zbasis, left_S_set_action, e_to_Se,
            extra_greedy=False, ensure_minimal=True, verbose=False)
        root = ResolutionNode(self, mod0, None, mat0)
        self.root = root
        self.node_cache = {}

    def make_actions(self, module):
        """Given a list of idempotents [e1, ..., en] representing ZSe1(+)...(+)ZSen,
        Return a list of "action" vectors, one for each element s of S, such that
        the Z-basis element i gets sent to the Z-basis element action[i]
        under multiplication by s.
        """
        e_to_Se = self.e_to_Se
        e_to_s_to_ii = self.e_to_s_to_ii
        op = self.op
        N = sum(map(len, map(e_to_Se.get, module)))
        offset = 0
        S_action_table = [[None] * N for _ in range(len(op))]
        for e in module:
            Se = e_to_Se[e]
            se_to_ii = e_to_s_to_ii[e]
            index_se_pairs = list(enumerate(Se, offset))
            for op_s1, table_s1 in zip(op, S_action_table):
                for index, se in index_se_pairs:
                    table_s1[index] = offset + se_to_ii[op_s1[se]]
            offset += len(Se)
        return list(map(Vector, S_action_table))

    def homology_list(self, maxdim, *, right_S_set_action=None, check=True, verbose=False):
        op = self.op
        S = range(len(op))
        if right_S_set_action is None:
            right_S_set_action = [[0] * len(op)]
        X = range(len(right_S_set_action))
        if check:
            act = right_S_set_action
            for row in right_S_set_action:
                if len(row) != len(S):
                    raise ValueError(f"len(right_S_set_action[i])={len(row)} mismatches {len(op)=}")
            for x in X:
                for i in S:
                    xi = act[x][i]
                    for j in S:
                        ij = op[i][j]
                        if act[xi][j] != act[x][ij]:
                            raise ValueError(f"act[act[{x}][{i}]][{j}]=act[{xi}][{j}]={act[xi][j]}, but"
                                             f"act[{x}][op[{i}][{j}]]=act[{x}][{ij}]={act[x][ij]}")
        e_to_Xe = {e: sorted(right_S_set_action[x][e] for x in X)
                   for e in self.e_to_Se}
        e_to_x_to_ii = {
            e: {x: ii for ii, x in enumerate(Xe)}
            for e, Xe in e_to_Xe.items()
        }

        @cache
        def outgoing_invariants(node) -> list:
            return node.outgoing_tensored_invariants(right_S_set_action, e_to_Xe, e_to_x_to_ii)

        @cache
        def homology(node) -> Counter:
            incoming = Counter()
            for child in node.get_children(verbose=verbose):
                incoming.update(outgoing_invariants(child))
            chains_rank = sum(map(len, map(e_to_Xe.get, node.module)))
            outgoing_rank = len(outgoing_invariants(node))
            incoming_rank = incoming.total()
            free_rank = chains_rank - outgoing_rank - incoming_rank
            result = Counter(d for d in incoming if d > 1) # torsion
            result[0] = free_rank
            return result

        @cache
        def homology_with_shift(node, shift) -> Counter:
            if shift == 0:
                return homology(node)
            else:
                result = Counter()
                for child in node.get_children(verbose=verbose):
                    result += homology_with_shift(child, shift - 1)
                return result

        return [invariant_factors(homology_with_shift(self.root, dim))
                for dim in range(maxdim + 1)]


class ResolutionNode:
    __slots__ = [
        "resolution", # The ProjectiveResolution we belong to
        "module", # a list of idempotents [e1, ..., en]: this node represents ZSe1 (+) ... (+) ZSen
        "prev_module", # The idempotents for the node one dimension lower; the target of the outgoing map
        "e_images",
        "children", # A list of nodes one dimension higher used to cover the kernel of the outgoing map
    ]
    def __init__(self, resolution: ProjectiveResolution, module, prev_module, e_images):
        self.resolution = resolution
        self.module = module
        self.prev_module = prev_module
        self.e_images = e_images
        self.children = None
        assert len(e_images) == len(module)

    def get_children(self, *, verbose=False, max_size_to_ensure_minimal=1000, max_size_to_cache=1000, max_size_for_extra_greedy=200):
        if self.children is not None:
            return self.children

        e_to_Se = self.resolution.e_to_Se
        if self.prev_module is None:
            assert self is self.resolution.root
            S_action_on_image = self.resolution.left_S_set_action
        else:
            S_action_on_image = self.resolution.make_actions(self.prev_module)
        # Forget the S-module structure and look at the Z-matrix
        outgoing_matrix_columns = [
            vec.shuffled_by_action(S_action_on_image[se])
            for vec, e in zip(self.e_images, self.module)
            for se in e_to_Se[e]
        ]
        kernel_basis = self.resolution.kernel_implementation(outgoing_matrix_columns)

        # Conversion between Z-basis indexes and summand ZSe indexes
        index_to_gen_index = []
        gen_index_to_index_range = []
        for i, e in enumerate(self.module):
            n0 = len(index_to_gen_index)
            index_to_gen_index.extend([i] * len(self.resolution.e_to_Se[e]))
            n1 = len(index_to_gen_index)
            gen_index_to_index_range.append(list(range(n0, n1)))

        # See if the kernel splits along the given summands
        indexes, summands = Lattice(len(outgoing_matrix_columns), kernel_basis).decompose(gen_index_to_index_range)
        if verbose:
            print(f"split into {len(indexes)} bins: {[len(x) for x in indexes]}")
        children = []
        for index_group, summand in zip(indexes, summands):
            gen_indexes = sorted({index_to_gen_index[ix] for ix in index_group})
            summand_gens = [self.module[gen_index] for gen_index in gen_indexes]
            split_e_images, split_next_module = cover_submodule_with_actions(
                summand.get_basis(),
                self.resolution.make_actions(summand_gens),
                e_to_Se,
                extra_greedy=(summand.rank <= max_size_for_extra_greedy),
                ensure_minimal=(summand.rank <= max_size_to_ensure_minimal),
                verbose=verbose,
            )
            if sum(map(len, split_e_images)) > max_size_to_cache:
                child = ResolutionNode(self.resolution, split_next_module, summand_gens, split_e_images)
            else:
                cache_key = tuple(split_next_module), tuple(summand_gens), tuple(map(tuple, split_e_images))
                node_cache = self.resolution.node_cache
                child = node_cache.get(cache_key)
                if child is None:
                    child = ResolutionNode(self.resolution, split_next_module, summand_gens, split_e_images)
                    node_cache[cache_key] = child
                else:
                    if verbose:
                        print(f"cache hit on {len(summand_gens)}gens <-- {len(split_next_module)}gens")
            children.append(child)
        self.children = children
        return children

    def outgoing_tensored_invariants(self, right_S_set_action, e_to_Xe, e_to_x_to_ii) -> list[int]:
        if self.prev_module is None:
            assert self is self.resolution.root
            # This is the root--the outgoing map is deleted.
            return []

        X = range(len(right_S_set_action))
        e_to_Se = self.resolution.e_to_Se
        # e_to_s_to_ii = self.resolution.e_to_s_to_ii

        # Tensoring moves from a rank-N0 module sum(ZSej)
        #                   to a rank-N1 module sum(ZXej)
        N0 = sum(map(len, map(e_to_Se.get, self.prev_module)))
        assert set(map(len, self.e_images)) <= {N0}
        N1 = sum(map(len, map(e_to_Xe.get, self.prev_module)))

        e_images = self.e_images

        def make_X_actions():
            # "actions" to move from sum(ZSe) to sum(ZXe).
            offset0 = offset1 = 0
            # the table that takes (index in X, index in union of Se) --> index in union of Xe
            x_action_table = [[None] * N0 for x in X]
            for e in self.prev_module:
                Se = e_to_Se[e]
                Xe = e_to_Xe[e]
                x_to_ii = e_to_x_to_ii[e]
                index0_se_pairs = list(enumerate(Se, offset0))
                for act_x, table_x in zip(right_S_set_action, x_action_table):
                    for index0, se in index0_se_pairs:
                        table_x[index0] = offset1 + x_to_ii[act_x[se]]
                offset0 += len(Se)
                offset1 += len(Xe)
            return list(map(Vector, x_action_table))
        x_actions = make_X_actions()

        def make_tensored_image():
            tensored_image = Lattice(N1, maxrank=len(X)*len(self.module))
            for e_image in e_images:
                for x_act in x_actions:
                    acted = e_image.shuffled_by_action(x_act, N1)
                    tensored_image.add_vector(acted)
            return tensored_image
        tensored_image = make_tensored_image()
        return tensored_image.nonzero_invariants()


