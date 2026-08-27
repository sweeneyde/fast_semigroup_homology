from collections import deque
from itertools import combinations
import string
from math import comb
import random

from .crs import CRS
from .knuth_bendix import kb_normalize

SYMBOLS = string.ascii_uppercase + string.ascii_lowercase + string.digits

def representation_by_generators(op, gens):
    """Given an operation table and a set of generators,
    represent all monoid elements using the generators.
    """
    n = len(op)
    rn = range(n)
    [identity] = [e for e in rn if all(op[e][x] == x == op[x][e] for x in rn)]
    representation = [None] * n
    q = deque([((), identity)])
    while q:
        tup, x = q.popleft()
        if representation[x] is None:
            representation[x] = tup
            for y in gens:
                q.append((tup + (y,), op[x][y]))
    return representation

def sample_permutations(arr, num_perms):
    if num_perms == 1:
        return [arr]
    elif num_perms == 2:
        return [arr, arr[::-1]]
    else:
        R = random.Random(0)
        arr = list(arr)
        result = set()
        for _ in range(num_perms):
            R.shuffle(arr)
            result.add(tuple(arr))
        return result

def crs_from_representation(op, gens, rep):
    alphabet = SYMBOLS[:len(gens)]
    gen_to_symbol = dict(zip(gens, alphabet))
    elements = [
        "".join(map(gen_to_symbol.get, rep[a]))
        for a in range(len(op))
    ]
    relations = [
         ("".join(map(gen_to_symbol.get, rep[a] + rep[b])),
          "".join(map(gen_to_symbol.get, rep[op[a][b]])))
        for a in range(len(op))
        for b in range(len(op))
    ]
    relations.sort()
    alphabet, relations = kb_normalize(alphabet, relations)
    crs = CRS(alphabet, relations)
    return elements, crs

def generate_crs_alternatives_with_generators(op, gens, perms_per_combo):
    orders = sample_permutations(gens, perms_per_combo)
    for gens1 in orders:
        rep = representation_by_generators(op, gens1)
        for gens2 in orders:
            elements, crs = crs_from_representation(op, gens2, rep)
            print(gens, gens1, gens2, "-->", crs.essential_counts(5))
            yield elements, crs

def generate_crs_alternatives(op, max_generators, max_combos, perms_per_combo):
    n = len(op)
    if max_generators is None:
        max_generators = n - 1
    rn = range(n)
    [identity] = [e for e in rn if all(op[e][x] == x == op[x][e] for x in rn)]
    non_identity = sorted(set(rn) - {identity})
    for k in range(max_generators + 1):
        if max_combos is None or comb(len(non_identity), k) <= max_combos:
            it = combinations(non_identity, k)
        else:
            R = random.Random(0)
            it = (R.sample(non_identity, k) for _ in range(max_combos))
        for gens in it:
            rep = representation_by_generators(op, gens)
            if None in rep:
                continue
            yield from generate_crs_alternatives_with_generators(op, gens, perms_per_combo)

def find_best_gens_crs(op, maxdim=5, *, max_generators=None, max_combos=None, perms_per_combo=1):
    def cost(elements_crs):
        _, crs = elements_crs
        return max(crs.essential_counts(maxdim))
    elements, crs = min(generate_crs_alternatives(op, max_generators, max_combos, perms_per_combo), key=cost)
    elements = list(map(crs.reduce, elements))
    return elements, crs
