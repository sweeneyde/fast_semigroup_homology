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
        root = EssentialCell.make_root()
        self._essentials = [(root,),
                            tuple(map(root.make_essential_1cell, self.alphabet))]
        self._coboundary_matrix_cache = {}
        self._cohomology_cache = {}

    def compute_essentials(self, maxdim):
        while len(self._essentials) <= maxdim:
            ess = []
            for x in self._essentials[-1]:
                ess.extend(x.make_essential_extensions(self))
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
            return (CollapsibleCell, None)
        for j in range(1, len(b)):
            if self.reducible(a + b[:j]):
                return (RedundantCell, (self.intern(b[:j]),
                                        self.intern(b[j:])))
        return (EssentialCell, None)

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
            for face in cell.get_faces(self)
            if type(face) is RedundantCell
        ]
        while stack:
            peek = stack[-1]
            peek._represent_inessrep()
            if peek._inessrep:
                stack.extend(peek._inessrep)
            else:
                stack.pop()

        result = []
        for cell in self._essentials[dim]:
            d = defaultdict(int)
            for face in cell._faces[0::2]:
                # even index --> don't negate
                if type(face) is EssentialCell:
                    d[face] += 1
                elif type(face) is RedundantCell:
                    assert not face._inessrep
                    for rep, q in face._essrep.items():
                        d[rep] += q
            for face in cell._faces[1::2]:
                # odd index --> do negate
                if type(face) is EssentialCell:
                    d[face] -= 1
                elif type(face) is RedundantCell:
                    assert not face._inessrep
                    for rep, q in face._essrep.items():
                        d[rep] -= q
            for k in [k for k, v in d.items() if not v]:
                del d[k]
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
            cell = cell.get_extension(x, self)
            if cell is None:
                return {}
        tp = type(cell)
        if tp is EssentialCell:
            return {tup: 1}
        elif tp is RedundantCell:
            essrep = cell.essential_representation(reduce=True)
            return {rep.to_tuple(): qty for rep, qty in essrep.items()}
        else:
            return {}

CRS = CompleteRewritingSystem

# Idea: don't represent the COLLAPSIBLE and DEGENERATE cells at all.
# We could even have two classes, on for ESSENTIAL and one for REDUNDANT.

extension_cache_hits = Counter()
extension_cache_misses = Counter()
face_cache_hits = Counter()
face_cache_misses = Counter()

class _CellBase:
    __slots__ = [
        "all_but_last", # Cell | None
        "last", # str | None
        "_extensions", # dict[str, Cell]
        "_faces", # None | tuple[Cell]
    ]

    def to_tuple(self):
        parts = []
        while (x := self.last) is not None:
            parts.append(x)
            self = self.all_but_last
        parts.reverse()
        return tuple(parts)

    def __str__(self):
        return "|".join(self.to_tuple()).join("[]")

    def __repr__(self):
        return f"{type(self).__name__}('{self!s}')"

    def get_faces(self, crs):
        if self._faces is not None:
            return self._faces
        last = self.last
        abl = self.all_but_last
        faces = (
            *[None if f is None else f.get_extension(last, crs)
              for f in abl.get_faces(crs)[:-1]],
            abl.all_but_last.get_extension(
                crs.operation(abl.last, last), crs
            ),
            abl,
        )
        self._faces = faces
        return faces

class EssentialCell(_CellBase):
    __slots__ = ()

    @staticmethod
    def make_root():
        root = EssentialCell()
        root.all_but_last = root.last = None
        root._extensions = {}
        root._faces = ()
        return root

    def make_essential_1cell(self, nxt):
        assert self.last is None
        child = EssentialCell()
        child.all_but_last = self
        child.last = nxt
        child._extensions = {}
        child._faces = (self, self)
        self._extensions[nxt] = child
        return child

    def make_essential_extension(self, nxt):
        child = EssentialCell()
        child.all_but_last = self
        child.last = nxt
        child._extensions = {}
        child._faces = None
        self._extensions[nxt] = child
        return child

    def make_essential_extensions(self, crs):
        return map(self.make_essential_extension,
                   crs._successor_words[self.last])

    def get_essential_extension(self, nxt):
        extensions = self._extensions
        if (cached := extensions.get(nxt)) is not None:
            return cached
        result = self.make_essential_extension(nxt)
        extensions[nxt] = result
        return result

    def get_extension(self, nxt, crs):
        extensions = self._extensions
        if (cached := extensions.get(nxt)) is not None:
            return cached
        if not nxt:
            return None
        _faces = None
        if self.last is None:
            _faces = (self, self)
            head = crs.intern(nxt[0])
            tail = crs.intern(nxt[1:])
            assert tail
            child = RedundantCell()
            child.paired_coface = self._extensions[head].get_extension(tail, crs)
        else:
            tp, splitting = crs.classify_pair(self.last, nxt)
            if tp is RedundantCell:
                pre, suf = splitting
                child = RedundantCell()
                child.paired_coface = self.get_essential_extension(pre).get_extension(suf, crs)
            else:
                child = tp()
        child.all_but_last = self
        child.last = nxt
        child._extensions = {}
        child._faces = _faces
        extensions[nxt] = child
        if type(child) is RedundantCell:
            child._set_inessrep_and_essrep(crs)
        return child

class CollapsibleCell(_CellBase):
    __slots__ = ()
    def get_extension(self, nxt, _):
        extensions = self._extensions
        if (cached := extensions.get(nxt)) is not None:
            return cached
        if not nxt:
            return None
        child = CollapsibleCell()
        child.all_but_last = self
        child.last = nxt
        child._extensions = {}
        child._faces = None
        extensions[nxt] = child
        return child

class RedundantCell(_CellBase):
    __slots__ = [
        "_essrep",
        "_inessrep",
        "paired_coface", # Cell
    ]

    def get_extension(self, nxt, crs):
        extensions = self._extensions
        if (cached := extensions.get(nxt)) is not None:
            return cached
        if not nxt:
            return None
        child = RedundantCell()
        child.all_but_last = self
        child.last = nxt
        child._extensions = {}
        child._faces = None
        child.paired_coface = self.paired_coface.get_extension(nxt, crs)
        extensions[nxt] = child
        child._set_inessrep_and_essrep(crs)
        return child

    def _set_inessrep_and_essrep(self, crs):
        faces = self.paired_coface.get_faces(crs)
        i = faces.index(self)
        bit = i & 1
        essrep = self._essrep = defaultdict(int)
        inessrep = self._inessrep = defaultdict(int)
        # Opposite parity face --> multiplicity 1
        for pos in faces[1-bit::2]:
            if type(pos) is EssentialCell:
                essrep[pos] += 1
            elif type(pos) is RedundantCell:
                inessrep[pos] += 1        
        # Same parity face --> multiplicty -1
        for neg in faces[bit:i:2] + faces[i+2::2]:
            if type(neg) is EssentialCell:
                essrep[neg] -= 1
            elif type(neg) is RedundantCell:
                inessrep[neg] -= 1

    def _represent_inessrep(self):
        if not self._inessrep:
            return
        old_inessrep = self._inessrep
        inessrep = self._inessrep = defaultdict(int)
        essrep = self._essrep
        for rep, q in old_inessrep.items():
            if not q:
                continue
            for reprep, qq in rep._essrep.items():
                essrep[reprep] += q*qq
            for reprep, qq in rep._inessrep.items():
                inessrep[reprep] += q*qq

    def essential_representation(self, reduce=False):
        if reduce:
            while self._inessrep:
                self._represent_inessrep()
        else:
            if self._inessrep:
                raise ValueError("still needs reducing")
        return dict(self._essrep)
