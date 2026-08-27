from collections import defaultdict

from mutable_lattice import Vector

from ..projective_resolution import ResolutionNode, ProjectiveResolution
from .homology_with_generators import homology

class ResolutionWithBarComparison:
    """
    For a finite monoid S, store the data of a projective resolution
    of the trivial ZS-module Z, along with the data of a chain map
    to the bar resolution.
    """
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
    """
    Represent a map from a ResolutionNode to the bar resolution.
    """
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
            # If the base_node ZSe is the root, map it to ZS[] in the bar resolution
            assert base_node.resolution.root is base_node
            assert len(base_node.module) == 1
            self.e_images_in_bar = [{(): 1}]
            assert index_as_child is None
            return
        assert base_node is prev.base_node.children[index_as_child]
        gen_indexes = prev.base_node.child_gen_indexes[index_as_child]
        assert [prev.base_node.module[e] for e in gen_indexes] == base_node.prev_module
        e_to_Se = base_node.resolution.e_to_Se
        identity = base_node.resolution.identity

        # If not at the root, define the chain map f on each generator x by
        #    f(e) := h(f(d(e))
        #       where d is the boundary (represented by base_node.e_images),
        #       f is already defined on the previous module,
        #       and h is defined by h(x0[x1|...|xn]) = 1[x0|x1|...|xn],
        # The map h satisfies dh+hd = id, so we can use it to lift a cycle
        # in the bar resoltution to one of its coboundaries: assuming d(x)=0,
        # we have d(h(x))=d(h(x))+h(d(x))=(dh+hd)(x)=0, so h(x) is a coboundary of x.
        # Because d(f(d(e))) = f(d(d(e))) = f(0) = 0, f(d(e)) is indeed a cycle,
        # so we can lift using h to get a valid chain map: d(f(e)) = d(h(f(d(e))) = f(d(e)).
        # The listed e_images_in_bar represents the vector 1[x1|...|xn] that
        # the 1 in ZS would get sent to if present. If only ZSe is present,
        # its e gets sent to e[x1|...|xn].
        e_images_in_bar = []
        for e_image in base_node.e_images:
            # This is a "very sparse vector", mapping tuple[int] --> int.
            bar_image = defaultdict(int)
            offset = 0
            # Iterate over the entries of d(e), one summand at a time
            for gi, e in zip(gen_indexes, base_node.prev_module, strict=True):
                prev_bar_image = prev.e_images_in_bar[gi]  # f(this part of d(e))
                Se = e_to_Se[e]
                for ii in range(len(Se)):
                    coeff1 = e_image[offset + ii]
                    if coeff1:
                        # For each nonzero entry of d(e), add in the corresponding
                        # multiple of the image of h(f(...))
                        for tup, coeff2 in prev_bar_image.items():
                            bar_image[(Se[ii],) + tup] += coeff1 * coeff2
                offset += len(Se)
            assert offset == len(e_image)
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
        h = homology(incoming, outgoing)
        generators = h.get_generators()
        invariants = h.get_invariants()
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
