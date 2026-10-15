"""
Implement the "collapsing scheme" defined by

    Brown, K.S. (1992).
    The Geometry of Rewriting Systems: A Proof of the Anick-Groves-Squier Theorem.
    In: Baumslag, G., Miller, C.F. (eds) Algorithms and Classification in Combinatorial Group Theory.
    Mathematical Sciences Research Institute Publications, vol 23. Springer, New York, NY.
    https://doi.org/10.1007/978-1-4613-9730-4_6

for the nerve of a monoid with a complete rewriting system.
"""

from collections import defaultdict, Counter, deque
from enum import Enum
from itertools import permutations

from mutable_lattice import Vector, Lattice, transpose

from ..homology_with_generators import homology

class CellKind(Enum):
    DEGENERATE = 0
    ESSENTIAL = 1
    REDUNDANT = 2
    COLLAPSIBLE = 3

DEGENERATE, ESSENTIAL, REDUNDANT, COLLAPSIBLE = CellKind

class CompleteRewritingSystem:
    __slots__ = [
        "alphabet",
        "rules",
        "max_rewrites",
        "_interned",
        "_operation_cache",
        "_successor_words",
        "_essentials",
        "_classify_pair_cache",
        "_coboundary_matrix_cache",
        "_cohomology_cache",
    ]

    def intern(self, word):
        return self._interned.setdefault(word, word)

    def reduce(self, word):
        word0 = word
        for _ in range(self.max_rewrites + 1):
            old = word
            for left, right in self.rules.items():
                word = word.replace(left, right)
            if word == old:
                return word
        raise RuntimeError(f"No fixed point was found for {word0!r} "
                           f"after {self.max_rewrites} loops. "
                           f"Stopped at {word!r}")

    def reducible(self, word):
        return any(left in word for left in self.rules.keys())

    def irreducible(self, word):
        return not self.reducible(word)

    def operation(self, a, b):
        ca = self._operation_cache[a]
        if (cached := ca.get(b)) is not None:
            return cached
        result = self.intern(self.reduce(a + b))
        ca[b] = result
        return result

    def __repr__(self):
        return f"CRS({self.alphabet!r}, {self.rules})"

    def __init__(self,
                 alphabet : str,
                 rules : dict[str, str],
                 max_rewrites=1000):
        alphabet = "".join(alphabet)
        rules = rules.copy()
        self.alphabet = alphabet
        self.rules = rules
        self.max_rewrites = max_rewrites

        for left, right in rules.items():
            if not isinstance(left, str) or not isinstance(right, str):
                raise ValueError("Rules must be dict[str,str]")
        set_alphabet = set(alphabet)
        if len(set_alphabet) != len(alphabet):
            raise ValueError("duplicate letters in alphabet")
        letters = {c for rule in rules.items() for side in rule for c in side}
        if not (letters <= set_alphabet):
            raise ValueError(f"{letters - set_alphabet} not in alphabet")
        for right in rules.values():
            for left in rules.keys():
                if left in right:
                    raise ValueError(f"left side {left!r} "
                                     f"contains right side {right!r}")
        for left1, left2 in permutations(rules.keys(), 2):
            if left1 in left2:
                raise ValueError(f"left side {left2!r} "
                                 f"contains left side {left1!r}")
        for left in rules.keys():
            if len(left) <= 1:
                raise ValueError(f"left side {left!r} is too small")

        self._interned = {}

        # Collect proper nonempty prefixes/suffixes
        pre_suf_to_right = defaultdict(dict)
        suf_pre_to_right = defaultdict(dict)
        for left, right in rules.items():
            for i in range(1, len(left)):
                pre, suf = left[:i], left[i:]
                pre_suf_to_right[pre][suf] = right
                suf_pre_to_right[suf][pre] = right
        # Look for critical pairs: compare (ab)c to a(bc).
        for b in pre_suf_to_right.keys() & suf_pre_to_right.keys():
            for c, bc in pre_suf_to_right[b].items():
                for a, ab in suf_pre_to_right[b].items():
                    if self.reduce(ab + c) != self.reduce(a + bc):
                        raise ValueError(f"Critical pair {(a, b, c)} did not resolve")

        stack = list(map(self.intern, alphabet))
        successor_words = {None: list(stack)}
        while stack:
            word = stack.pop()
            if word in successor_words:
                continue
            result = []
            for i in range(len(word)):
                for suf, right in pre_suf_to_right[word[i:]].items():
                    if self.irreducible((word + suf)[:-1]):
                        # word + suf reduces to right,
                        # but no proper prefix reduces.
                        result.append(self.intern(suf))
            successor_words[word] = tuple(result)
            stack.extend(result)

        self._operation_cache = defaultdict(dict)
        self._successor_words = successor_words
        self._classify_pair_cache = {}
        self._essentials = [(Cell.make_root(self),)]
        self._coboundary_matrix_cache = {}
        self._cohomology_cache = {}

    def compute_essentials(self, maxdim):
        while len(self._essentials) <= maxdim:
            ess = []
            for x in self._essentials[-1]:
                ess.extend(x.get_essential_extensions())
            self._essentials.append(tuple(ess))

    def essential_tuples(self, maxdim):
        self.compute_essentials(maxdim)
        to_tuple = {self._essentials[0][0]: ()}
        for ess_list in self._essentials[1:maxdim+1]:
            for ess in ess_list:
                to_tuple[ess] = to_tuple[ess.all_but_last] + (ess.last,)
        return [list(map(to_tuple.get, ess_list))
                for ess_list in self._essentials[:maxdim+1]]

    def essential_counts(self, maxdim):
        assert maxdim >= 1
        counts = [1, len(self.alphabet)]
        last_counts = Counter(self.alphabet)
        while len(counts) <= maxdim:
            next_last_counts = Counter()
            for last, count in last_counts.items():
                for word in self._successor_words[last]:
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
                return (REDUNDANT, (self.intern(b[:j]),
                                    self.intern(b[j:])))
        return (ESSENTIAL, None)

    def classify_pair(self, a, b):
        if (cached := self._classify_pair_cache.get((a, b))) is not None:
            return cached
        result = self._classify_pair_internal(a, b)
        self._classify_pair_cache[a, b] = result
        return result

    def boundary_dicts(self, dim):
        if dim == 0:
            return [{}]
        self.compute_essentials(dim)

        # Calculate essential representations of all faces
        # of all essential cells
        stack = [
            face
            for cell in self._essentials[dim]
            for face in cell.get_faces()
            if face.kind is REDUNDANT
        ]
        while stack:
            cell = stack.pop()
            if cell._inessrep:
                cell._update_essrep()
                if cell._inessrep:
                    stack.append(cell)
                    stack.extend(cell._inessrep.keys())

        result = []
        for cell in self._essentials[dim]:
            d = defaultdict(int)
            for i, face in enumerate(cell._faces):
                sign = (-1)**i
                if face.kind is ESSENTIAL:
                    d[face] += sign
                elif face.kind is REDUNDANT:
                    for rep, qty in face._essrep.items():
                        d[rep] += sign*qty
            result.append(dict(d))
        return result

    def boundary_vectors(self, dim):
        if dim == 0:
            return [Vector([])]
        self.compute_essentials(dim)
        lower_essentials = self._essentials[dim - 1]
        N = len(lower_essentials)
        cell_to_index = {cell: i for i, cell in enumerate(lower_essentials)}
        result = []
        for d in self.boundary_dicts(dim):
            v = Vector.zero(N)
            for rep, qty in d.items():
                v[cell_to_index[rep]] = qty
            result.append(v)
        return result

    def coboundary_matrix(self, dim):
        if dim == -1:
            return []
        if dim not in self._coboundary_matrix_cache:
            boundary = self.boundary_vectors(dim + 1)
            N = len(self._essentials[dim])
            self._coboundary_matrix_cache[dim] = transpose(N, boundary)
        return self._coboundary_matrix_cache[dim]

    def boundary_nonzero_invariants(self, dim):
        if dim == 0:
            return []
        N = len(self._essentials[dim - 1])
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
            free_rank = len(self._essentials[dim]) - len(outgoing) - len(incoming)
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
            free_rank = len(self._essentials[dim]) - len(outgoing) - len(incoming)
            torsion = dict(Counter(outgoing))
            torsion.pop(1, None)
            if free_rank:
                result.append({0 : free_rank} | torsion)
            else:
                result.append(torsion)
        return result

    def cohomology_with_generators(self, dim):
        if dim not in self._cohomology_cache:
            incoming = self.coboundary_matrix(dim - 1)
            outgoing = self.coboundary_matrix(dim)
            h = homology(incoming, outgoing)
            generators = h.get_generators()
            invariants = h.get_invariants()
            ess_dim = self._essentials[dim]
            bar_gens = [{ess_dim[i].to_tuple() : gen[i]
                         for i in filter(gen.__getitem__, range(len(ess_dim)))}
                        for gen in generators]
            self._cohomology_cache[dim] = invariants, bar_gens
        return self._cohomology_cache[dim]

    def essential_representation(self, tup):
        [cell] = self._essentials[0]
        for x in tup:
            cell = cell.get_extension(x)
        if cell.kind is ESSENTIAL:
            return {tup: 1}
        elif cell.kind is REDUNDANT:
            while cell._inessrep:
                cell._update_essrep()
            return {rep.to_tuple(): qty for rep, qty in cell._essrep.items()}
        else:
            return {}

CRS = CompleteRewritingSystem

# Idea: don't represent the COLLAPSIBLE and DEGENERATE cells at all.
# We could even have two classes, on for ESSENTIAL and one for REDUNDANT.

class Cell:
    __slots__ = [
        "crs",
        "kind", # CellKind
        "all_but_last", # Cell | None
        "last", # str | None
        "_extensions", # dict[str, Cell]
        "_essential_extensions", # tuple[Cell]
        "_faces", # tuple[Cell] | None
        "paired_coface", # None | Cell
        "_inessrep", # dict[Cell, int]
        "_essrep", # dict[Cell, int]
    ]

    def to_tuple(self):
        parts = []
        while self.last is not None:
            parts.append(self.last)
            self = self.all_but_last
        parts.reverse()
        return tuple(parts)

    def __str__(self):
        return "|".join(reversed(self.to_tuple())).join("[]")

    def __repr__(self):
        return f"Cell({self.kind}, ..., {self.last})"

    @classmethod
    def make_root(cls, crs):
        root = cls()
        root.crs = crs
        root.kind = ESSENTIAL
        root.all_but_last = root.last = None
        root._extensions = {}
        root.paired_coface = None
        root._essential_extensions = None
        root._faces = ()
        root._inessrep = None
        root._essrep = None
        return root

    def make_child(self, last, kind, paired_coface=None):
        child = __class__()
        child.all_but_last = self
        child.last = last
        child.kind = kind
        child.paired_coface = paired_coface
        child.crs = self.crs
        child._extensions = {}
        child._faces = (self, self) if self.last is None else None
        child._essential_extensions = None
        self._extensions[last] = child
        if kind is REDUNDANT:
            essrep = defaultdict(int)
            inessrep = defaultdict(int)
            faces = paired_coface.get_faces()
            i = faces.index(child)
            for j, face in enumerate(faces):
                if i == j:
                    continue
                sign = (-1)**(i+j+1)
                if face.kind is ESSENTIAL:
                    essrep[face] += sign
                elif face.kind is REDUNDANT:
                    inessrep[face] += sign
            child._inessrep = inessrep
            child._essrep = essrep
        else:
            child._inessrep = None
            child._essrep = None
        return child

    def get_extension(self, nxt):
        extensions = self._extensions
        if (cached := extensions.get(nxt)) is not None:
            return cached
        mykind = self.kind
        paired_coface = None
        if self.last is None:
            assert mykind is ESSENTIAL
            if not nxt:
                kind = DEGENERATE
            elif len(nxt) == 1:
                kind = ESSENTIAL
            else:
                kind = REDUNDANT
                head = self.crs.intern(nxt[0])
                tail = self.crs.intern(nxt[1:])
                paired_coface = self.get_extension(head).get_extension(tail)
        else:
            kind = mykind
            if mykind is ESSENTIAL:
                kind, splitting = self.crs.classify_pair(self.last, nxt)
                if kind is REDUNDANT:
                    pre, suf = splitting
                    paired_coface = self.get_extension(pre).get_extension(suf)
            elif mykind is REDUNDANT:
                paired_coface = self.paired_coface.get_extension(nxt)
        new_cell = self.make_child(nxt, kind, paired_coface)
        return new_cell

    def get_essential_extensions(self):
        if (_ess := self._essential_extensions) is not None:
            return _ess
        if self.kind is not ESSENTIAL:
            raise ValueError
        succ = self.crs._successor_words[self.last]
        _ess = [self.get_extension(nxt) for nxt in succ]
        self._essential_extensions = _ess
        return _ess

    def get_faces(self):
        if self._faces is not None:
            return self._faces
        last = self.last
        abl = self.all_but_last
        faces = [f.get_extension(last) for f in abl.get_faces()[:-1]]
        faces.append(
            abl.all_but_last.get_extension(
                self.crs.operation(abl.last, last)
            )
        )
        faces.append(abl)
        # faces = tuple(faces)
        self._faces = faces
        return faces

    def _update_essrep(self):
        essrep = self._essrep
        old_inessrep = self._inessrep
        # Just asserting we're not re-entrant
        self._essrep = self._inessrep = object()
        inessrep = defaultdict(int)
        for rep, rep_qty in old_inessrep.items():
            for reprep, reprep_qty in rep._essrep.items():
                essrep[reprep] += rep_qty * reprep_qty
            for reprep, reprep_qty in rep._inessrep.items():
                inessrep[reprep] += rep_qty * reprep_qty
        self._essrep = essrep
        self._inessrep = inessrep
