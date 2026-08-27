"""
Implement the "collapsing scheme" defined by

    Brown, K.S. (1992).
    The Geometry of Rewriting Systems: A Proof of the Anick-Groves-Squier Theorem.
    In: Baumslag, G., Miller, C.F. (eds) Algorithms and Classification in Combinatorial Group Theory.
    Mathematical Sciences Research Institute Publications, vol 23. Springer, New York, NY.
    https://doi.org/10.1007/978-1-4613-9730-4_6

for the nerve of a monoid with a complete rewriting system.
"""


from collections import Counter, defaultdict, deque
from itertools import permutations
from enum import Enum

from mutable_lattice import Vector, Lattice, transpose
from ..homology_with_generators import homology

class CellType(Enum):
    DEGENERATE = 0
    ESSENTIAL = 1
    REDUNDANT = 2
    COLLAPSIBLE = 3

DEGENERATE = CellType.DEGENERATE
ESSENTIAL = CellType.ESSENTIAL
REDUNDANT = CellType.REDUNDANT
COLLAPSIBLE = CellType.COLLAPSIBLE

class CompleteRewritingSystem:

    __slots__ = [
        "alphabet",
        "rules",
        "max_rewrites",
        "essentials",
        "successor_words",
        "operation_cache",
        "classify_pair_cache",
        "essential_representation_cache",
        "coboundary_matrix_cache",
        "cohomology_cache",
    ]

    def reduce(self, word):
        word0 = word
        for _ in range(self.max_rewrites + 1):
            old = word
            for left, right in self.rules:
                word = word.replace(left, right)
            if word == old:
                return word
        raise RuntimeError(f"No fixed point was found for {word0!r} "
                           f"after {self.max_rewrites} loops. "
                           f"Stopped at {word!r}")

    def reducible(self, word):
        return any(left in word for left, _ in self.rules)

    def irreducible(self, word):
        return not self.reducible(word)

    def operation(self, a, b):
        if (cached := self.operation_cache.get((a, b))) is not None:
            return cached
        result = self.reduce(a + b)
        self.operation_cache[a, b] = result
        return result

    def __repr__(self):
        return f"CRS({self.alphabet!r}, {self.rules})"

    def __init__(self, alphabet, rules, max_rewrites=1000):
        if len(set(alphabet)) != len(alphabet):
            raise ValueError("alphabet has duplicates")
        self.alphabet = "".join(alphabet)
        self.rules = tuple(rules)
        self.max_rewrites = max_rewrites
        # Assert we're provided correct arguments
        for left, right in rules:
            if not isinstance(left, str) or not isinstance(right, str):
                raise ValueError("Rules must be a list of pairs of strings")
        set_alphabet = set(alphabet)
        if len(set_alphabet) != len(alphabet):
            raise ValueError("duplicate letters in alphabet")
        letters = {c for rule in rules for side in rule for c in side}
        if not (letters <= set_alphabet):
            raise ValueError(f"{letters - set_alphabet} not in alphabet")
        for _, right in rules:
            for left, _ in rules:
                if left in right:
                    raise ValueError(f"left side {left!r} "
                                     f"contains right side {right!r}")
        for (left1, _), (left2, _) in permutations(rules, 2):
            if left1 in left2:
                raise ValueError(f"left side {left2!r} "
                                 f"contains left side {left1!r}")
        for left, _ in rules:
            if len(left) <= 1:
                raise ValueError(f"left side {left!r} is too small")
        # Collect proper nonempty prefixes/suffixes
        prefix_to_rules = defaultdict(list)
        suffix_to_rules = defaultdict(list)
        for left, right in rules:
            for i in range(1, len(left)):
                prefix, suffix = left[:i], left[i:]
                prefix_to_rules[prefix].append((suffix, right))
                suffix_to_rules[suffix].append((prefix, right))
        # Look for critical pairs: comapare (ab)c to a(bc).
        for middle in prefix_to_rules.keys() & suffix_to_rules.keys():
            prefix_rules = prefix_to_rules[middle]
            suffix_rules = suffix_to_rules[middle]
            for suffix_rule_prefix, suffix_rule_right in suffix_rules:
                for prefix_rule_suffix, prefix_rule_right in prefix_rules:
                    result1 = self.reduce(suffix_rule_right + prefix_rule_suffix)
                    result2 = self.reduce(suffix_rule_prefix + prefix_rule_right)
                    if result1 != result2:
                        raise ValueError(f"Critical pair {(suffix_rule_prefix, middle, prefix_rule_suffix)} did not resolve")

        self.essentials = [((),),
                           tuple([(a,) for a in alphabet])]

        stack = list(alphabet)
        successor_words = {}
        while stack:
            word = stack.pop()
            if word in successor_words:
                continue
            result = []
            for i in range(len(word)):
                suffix = word[i:]
                for left_suffix, right in prefix_to_rules.get(suffix, ()):
                    if self.irreducible((word + left_suffix)[:-1]):
                        # word + left_suffix reduces to right,
                        # but no proper prefix reduces.
                        result.append(left_suffix)
            successor_words[word] = tuple(result)
            stack.extend(result)

        self.successor_words = successor_words
        self.operation_cache = {}
        self.classify_pair_cache = {}
        self.essential_representation_cache = {}
        self.coboundary_matrix_cache = {}
        self.cohomology_cache = {}

    def compute_essentials(self, maxdim):
        while len(self.essentials) <= maxdim:
            ess = []
            starts = self.essentials[-1]
            for start in starts:
                for word in self.successor_words[start[-1]]:
                    ess.append(start + (word,))
            self.essentials.append(tuple(ess))

    def essential_counts(self, maxdim=0):
        counts = [1, len(self.alphabet)]
        last_counts = Counter(self.alphabet)
        while len(counts) <= maxdim:
            next_last_counts = Counter()
            for last, count in last_counts.items():
                for word in self.successor_words[last]:
                    next_last_counts[word] += count
            counts.append(next_last_counts.total())
            last_counts = next_last_counts
        return counts

    def elements(self, limit=1000):
        alphabet = self.alphabet
        q = deque([''])
        done = set()
        while q:
            x = q.popleft()
            if x in done:
                continue
            yield x
            done.add(x)
            if len(done) > limit:
                raise ValueError(f"found more than {limit=} elements")
            for a in alphabet:
                q.append(self.reduce(x + a))

    def op_table(self, elements):
        element_to_index = {x: i for i, x in enumerate(elements)}
        op = [[element_to_index[self.operation(x, y)] for y in elements]
              for x in elements]
        return op

    def _classify_pair_internal(self, a, b):
        # Evaluate a pair of adjacent entries of a cell.
        # COLLAPSIBLE cells have an adjacent pair concatenate to an irreducible.
        # REDUNDANT cells have an adjacent pair concatenate to have a reducible proper prefix.
        # ESSENTIAL cells have none of these.
        if self.irreducible(a + b):
            return (COLLAPSIBLE, None)
        for j in range(1, len(b)):
            if self.reducible(a + b[:j]):
                return (REDUNDANT, j)
        return (ESSENTIAL, None)

    def classify_pair(self, a, b):
        if (cached := self.classify_pair_cache.get((a, b))) is not None:
            return cached
        result = self._classify_pair_internal(a, b)
        self.classify_pair_cache[a, b] = result
        return result

    def classify_cell(self, cell):
        if len(cell) == 0:
            return (ESSENTIAL, None)
        if "" in cell:
            return (DEGENERATE, None)
        if len(cell[0]) > 1:
            collapsible = (cell[0][0], cell[0][1:]) + cell[1:]
            return (REDUNDANT, (collapsible, 1))
        for i in range(1, len(cell)):
            a, b = cell[i-1], cell[i]
            kind, index = self.classify_pair(a, b)
            if kind == ESSENTIAL:
                continue
            elif kind == COLLAPSIBLE:
                redundant = cell[:i-1] + (a + b,) + cell[i+1:]
                return (COLLAPSIBLE, (redundant, i))
            else:
                assert kind == REDUNDANT
                collapsible = cell[:i] + (b[:index], b[index:]) + cell[i+1:]
                return (REDUNDANT, (collapsible, i+1))
        return (ESSENTIAL, None)

    def get_faces(self, cell):
        for i in range(len(cell) + 1):
            if i == 0:
                face = cell[1:]
            elif i == len(cell):
                face = cell[:-1]
            else:
                x = self.operation(cell[i-1], cell[i])
                face = cell[:i-1] + (x,) + cell[i+1:]
            yield i, face

    def _essential_representation_internal(self, cell):
        kind, data = self.classify_cell(cell)
        if kind == ESSENTIAL:
            return {cell: 1}
        elif kind == COLLAPSIBLE:
            return {}
        elif kind == DEGENERATE:
            return {}
        assert kind == REDUNDANT
        collapsible, index = data
        result = defaultdict(int)
        for i, face in self.get_faces(collapsible):
            if i == index:
                continue
            sign = (-1)**(i +index + 1)
            for rep, count in self.essential_representation(face).items():
                result[rep] += sign*count
        return result

    def essential_representation(self, cell):
        if (cached := self.essential_representation_cache.get(cell)) is not None:
            return cached
        result = self._essential_representation_internal(cell)
        self.essential_representation_cache[cell] = result
        return result

    def essential_boundary(self, cell):
        assert self.classify_cell(cell) == (ESSENTIAL, None)
        result = defaultdict(int)
        for i, face in self.get_faces(cell):
            sign = (-1)**i
            for rep, count in self.essential_representation(face).items():
                result[rep] += sign*count
        return dict(result)

    def boundary_vectors(self, dim):
        if dim == 0:
            return []
        self.compute_essentials(dim)
        lower_essentials = self.essentials[dim - 1]
        N = len(lower_essentials)
        cell_to_index = {cell: i for i, cell in enumerate(lower_essentials)}
        result = []
        for cell in self.essentials[dim]:
            v = Vector.zero(N)
            for rep, count in self.essential_boundary(cell).items():
                v[cell_to_index[rep]] = count
            result.append(v)
        return result

    def coboundary_matrix(self, dim):
        if dim == -1:
            return []
        if dim not in self.coboundary_matrix_cache:
            boundary = self.boundary_vectors(dim + 1)
            N = len(self.essentials[dim])
            self.coboundary_matrix_cache[dim] = transpose(N, boundary)
        return self.coboundary_matrix_cache[dim]

    def boundary_nonzero_invariants(self, dim):
        if dim == 0:
            return []
        N = len(self.essentials[dim - 1])
        boundary_vectors = self.boundary_vectors(dim)
        L = Lattice(N, boundary_vectors, maxrank=len(boundary_vectors))
        return L.nonzero_invariants()

    def homology_list(self, maxdim):
        assert maxdim >= 0
        invariants = list(map(self.boundary_nonzero_invariants, range(maxdim + 2)))
        result = []
        for dim in range(maxdim + 1):
            outgoing = invariants[dim]
            incoming = invariants[dim + 1]
            free_rank = len(self.essentials[dim]) - len(outgoing) - len(incoming)
            torsion = dict(Counter(incoming))
            torsion.pop(1, None)
            if free_rank:
                result.append({0 : free_rank} | torsion)
            else:
                result.append(torsion)
        return result

    def cohomology_list(self, maxdim):
        invariants = list(map(self.boundary_nonzero_invariants, range(maxdim + 2)))
        result = []
        for dim in range(maxdim + 1):
            outgoing = invariants[dim]
            incoming = invariants[dim + 1]
            free_rank = len(self.essentials[dim]) - len(outgoing) - len(incoming)
            torsion = dict(Counter(outgoing))
            torsion.pop(1, None)
            if free_rank:
                result.append({0 : free_rank} | torsion)
            else:
                result.append(torsion)
        return result

    def cohomology_with_generators(self, dim):
        if dim not in self.cohomology_cache:
            incoming = self.coboundary_matrix(dim - 1)
            outgoing = self.coboundary_matrix(dim)
            h = homology(incoming, outgoing)
            generators = h.get_generators()
            invariants = h.get_invariants()
            ess_dim = self.essentials[dim]
            bar_gens = [{ess_dim[i] : gen[i]
                         for i in filter(gen.__getitem__, range(len(ess_dim)))}
                        for gen in generators]
            self.cohomology_cache[dim] = invariants, bar_gens
        return self.cohomology_cache[dim]

CRS = CompleteRewritingSystem
