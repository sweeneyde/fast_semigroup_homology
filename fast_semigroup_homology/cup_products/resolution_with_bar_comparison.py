from collections import defaultdict

from mutable_lattice import Vector

from ..projective_resolution import ResolutionNode, ProjectiveResolution
from .homology_with_generators import homology_with_generators

class ResolutionWithBarComparison:
    __slots__ = [
        "base_res", # ProjectiveResolution
        "nodes_by_dimension", # list[list[ComparisonNode]]
    ]
    def __init__(self, op, base_res=None):
        if base_res is None:
            base_res = ProjectiveResolution(op)
        else:
            if len(base_res.left_S_set_action[0]) != 1:
                raise ValueError("Expected a resolution of the trivial module Z")
        self.base_res = base_res
        self.nodes_by_dimension = [[BarComparisonNode(base_res.root, None, None)]]

    def extend_to_dimension(self, maxdim):
        nbd = self.nodes_by_dimension
        self.base_res.extend_to_dimension(maxdim)
        while len(nbd) - 1 < maxdim:
            next_dimension = []
            for comp_node in nbd[-1]:
                for i, child in enumerate(comp_node.base_node.children):
                    next_dimension.append(BarComparisonNode(child, comp_node, i))
            nbd.append(next_dimension)

    def homology_freepart_generators_in_bar(self, dim):
        self.extend_to_dimension(dim)
        results = []
        for comp_node in self.nodes_by_dimension[dim]:
            invariants, bar_generators = comp_node.homology_with_generators_in_bar()
            for d, bar_gen in zip(invariants, bar_generators):
                if d == 0:
                    results.append(bar_gen)
        return results

class BarComparisonNode:
    __slots__ = [
        # A ResolutionNode we're mapping out of
        "base_node",
        # Each summand's identity gets mapped to some 1[x1|...|xn] in the bar res.
        # We store a list[tuple[int]], with entries (x1,...,xn).
        # If the generator e is not the identity, then technically we map to e[x1|...|xn],
        # but we're only going to store the (x1,...,xn).
        "e_images_in_bar",
        # The BarComparisonNode our boundary goes to
        "prev",
        # self.base_node is self.parent.base_node.children[index_as_child]
        "index_as_child",
    ]

    def __init__(self, base_node : ResolutionNode, prev, index_as_child):
        self.base_node = base_node
        self.prev = prev
        self.index_as_child = index_as_child
        if prev is None:
            assert base_node.resolution.root is base_node
            assert len(base_node.module) == 1
            self.e_images_in_bar = [{(): 1}]
            assert index_as_child is None
            return
        assert base_node is prev.base_node.children[index_as_child]
        gen_indexes = prev.base_node.child_gen_indexes[index_as_child]
        assert [prev.base_node.module[e] for e in gen_indexes] == base_node.prev_module
        e_to_Se = base_node.resolution.e_to_Se

        e_images_in_bar = []
        for e_image in base_node.e_images:
            bar_image = defaultdict(int)
            offset = 0
            for gi, e in zip(gen_indexes, base_node.prev_module, strict=True):
                prev_bar_image = prev.e_images_in_bar[gi]
                Se = e_to_Se[e]
                for ii in range(len(Se)):
                    coeff1 = e_image[offset + ii]
                    if coeff1:
                        for tup, coeff2 in prev_bar_image.items():
                            bar_image[(Se[ii],) + tup] += coeff1 * coeff2
                offset += len(Se)
            assert offset == len(e_image)
            identity = base_node.resolution.identity
            e_images_in_bar.append({
                tup : count for tup, count in bar_image.items()
                if count
                if identity not in tup  # normalized bar resolution: no tuples containing id
            })

        self.e_images_in_bar = e_images_in_bar

    def homology_with_generators_in_bar(self):
        # To tensor a map, identify together each ZSe summand together to a single Z.
        e_to_Se = self.base_node.resolution.e_to_Se
        if self.prev is None:
            # Delete the augmentation map while tensoring
            outgoing = [Vector([]) for _ in self.base_node.e_images]
        else:
            outgoing_tensor_action = []
            for i, e in enumerate(self.base_node.prev_module):
                outgoing_tensor_action.extend([i] * len(e_to_Se[e]))
            outgoing_tensor_action = Vector(outgoing_tensor_action)
            outgoing = [e_image.shuffled_by_action(outgoing_tensor_action, len(self.base_node.prev_module))
                        for e_image in self.base_node.e_images]
        incoming = []
        for child, gen_indexes in zip(self.base_node.get_children(),
                                      self.base_node.child_gen_indexes, strict=True):
            incoming_tensor_action = []
            for gi, e in zip(gen_indexes, child.prev_module, strict=True):
                incoming_tensor_action.extend([gi] * len(e_to_Se[e]))
            incoming_tensor_action = Vector(incoming_tensor_action)
            for e_image in child.e_images:
                incoming.append(e_image.shuffled_by_action(incoming_tensor_action, len(outgoing)))
        invariants, generators = homology_with_generators(incoming, outgoing)
        # Push forward to the bar resolution
        bar_generators = []
        for gen in generators:
            bar_gen = defaultdict(int)
            for bar_image, coeff1 in zip(self.e_images_in_bar, gen):
                if not coeff1:
                    continue
                for tup, coeff2 in bar_image.items():
                    bar_gen[tup] += coeff1 * coeff2
            bar_generators.append({
                tup: count for tup, count in bar_gen.items()
                if count
            })
        return invariants, bar_generators
