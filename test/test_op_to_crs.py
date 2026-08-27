from fast_semigroup_homology.cup_products.rewriting.op_to_crs import (
    find_best_gens_crs,
)
import pytest

VERIFY_RULES = True

def check_crs(op, **kwargs):
    elements, crs = find_best_gens_crs(op, **kwargs)
    assert sorted(crs.elements()) == sorted(elements)
    assert set(crs.alphabet) <= set(elements)
    assert op == crs.op_table(elements)
    return elements, crs

@pytest.mark.parametrize("n", [2,3,4,5,6,7])
def test_cyclic(n):
    op = [[(i+j)%n for j in range(n)] for i in range(n)]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "A"
        assert elements == ["A"*i for i in range(n)]
        assert sorted(crs.rules) == [("A"*n, "")]
    assert crs.essential_counts(5) == [1, 1, 1, 1, 1, 1]

def test_V4():
    elements, crs = check_crs([[i^j for j in range(4)] for i in range(4)])
    if VERIFY_RULES:
        assert crs.alphabet == "AB"
        assert elements == ["", "A", "B", "AB"]
        assert sorted(crs.rules) == [("AA", ""), ("BA", "AB"), ("BB", "")]
    assert crs.essential_counts(5) == [1, 2, 3, 4, 5, 6]

def test_D6():
    op = [[0,1,2,3,4,5],
          [1,0,3,2,5,4],
          [2,5,4,1,0,3],
          [3,4,5,0,1,2],
          [4,3,0,5,2,1],
          [5,2,1,4,3,0]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "AB"
        assert elements == ["", "A", "AB", "B", "BA", "ABA"]
        assert sorted(crs.rules) == [("AA", ""), ("BAB", "ABA"), ("BB", "")]
    assert crs.essential_counts(5) == [1, 2, 3, 5, 9, 17]

def test_rect22():
    op = [[0,1,0,1,0],
          [0,1,0,1,1],
          [2,3,2,3,2],
          [2,3,2,3,3],
          [0,1,2,3,4]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "AB"
        assert elements == ["A", "AB", "BA", "B", ""]
        assert sorted(crs.rules) == [("AA", "A"), ("ABA", "A"), ("BAB", "B"), ("BB", "B")]
    assert crs.essential_counts(5) == [1, 2, 4, 8, 16, 32]

def test_0ZZZZ():
    op = [[0,1,0,1,0,0],
          [0,1,0,1,0,1],
          [2,3,2,3,2,2],
          [2,3,2,3,2,3],
          [0,1,0,1,0,4],
          [0,1,2,3,4,5]]

    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "ABC"
        assert elements == ["AC", "A", "BC", "B", "C", ""]
        assert sorted(crs.rules) == [('AA', 'A'), ('AB', 'A'), ('BA', 'B'), ('BB', 'B'), ('CA', 'A'), ('CB', 'A'), ('CC', 'AC')]
    assert crs.essential_counts(5) == [1, 3, 7, 15, 31, 63]

def test_rect22_adjoin_C2_trivially():
    op = [[0,1,0,1,0,0],
          [0,1,0,1,1,1],
          [2,3,2,3,2,2],
          [2,3,2,3,3,3],
          [0,1,2,3,4,5],
          [0,1,2,3,5,4]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "ABC"
        assert elements == ["A", "AB", "BA", "B", "", "C"]
        assert sorted(crs.rules) == [('AA', 'A'), ('ABA', 'A'), ('AC', 'A'), ('BAB', 'B'),
                                    ('BB', 'B'), ('BC', 'B'), ('CA', 'A'), ('CB', 'B'),
                                    ('CC', '')]
    assert crs.essential_counts(5) == [1, 3, 9, 27, 81, 243]

def test_rect22_adjoin_C2_right_action():
    op = [[0,1,0,1,0,1],
          [0,1,0,1,1,0],
          [2,3,2,3,2,3],
          [2,3,2,3,3,2],
          [0,1,2,3,4,5],
          [0,1,2,3,5,4]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "ABC"
        assert elements == ["A", "AC", "B", "BC", "", "C"]
        assert sorted(crs.rules) == [('AA', 'A'), ('AB', 'A'),
                                    ('BA', 'B'), ('BB', 'B'),
                                    ('CA', 'A'), ('CB', 'B'), ('CC', '')]
    assert crs.essential_counts(5) == [1, 3, 7, 15, 31, 63]

def test_rect22_adjoin_C2_left_action():
    op = [[0,1,0,1,0,0],
          [0,1,0,1,1,1],
          [2,3,2,3,2,2],
          [2,3,2,3,3,3],
          [0,1,2,3,4,5],
          [2,3,0,1,5,4]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "ABC"
        assert elements == ["A", "B", "CA", "CB", "", "C"]
        assert sorted(crs.rules) == [('AA', 'A'), ('AB', 'B'), ('AC', 'A'), ('BA', 'A'), ('BB', 'B'), ('BC', 'B'), ('CC', '')]
    assert crs.essential_counts(5) == [1, 3, 7, 15, 31, 63]

def test_rect22_adjoin_C2_2sided_action():
    op = [[0,1,0,1,0,1],
          [0,1,0,1,1,0],
          [2,3,2,3,2,3],
          [2,3,2,3,3,2],
          [0,1,2,3,4,5],
          [2,3,0,1,5,4]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "AB"
        assert elements == ["A", "AB", "BA", "BAB", "", "B"]
        assert sorted(crs.rules) == [('AA', 'A'), ('ABA', 'A'), ('BB', '')]
    assert crs.essential_counts(5) == [1, 2, 3, 5, 9, 17]


def test_rect22_adjoin_idempotent_duplicate():
    op = [[0,1,0,1,0,0],
          [0,1,0,1,0,1],
          [2,3,2,3,2,2],
          [2,3,2,3,2,3],
          [0,1,0,1,4,4],
          [0,1,2,3,4,5]]
    elements, crs = check_crs(op)
    if VERIFY_RULES:
        assert crs.alphabet == "AB"
        assert elements == ["BAB", "BA", "AB", "A", "B", ""]
        assert sorted(crs.rules) == [('AA', 'A'), ('ABA', 'A'), ('BB', 'B')]
    assert crs.essential_counts(5) == [1, 2, 3, 5, 9, 17]

def test_rect22_rect22():
    rect22 = [[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]]
    pairs = [(i, j) for i in range(5) for j in range(5)]
    pair_to_index = {ij: index for index, ij in enumerate(pairs)}
    op = [[pair_to_index[rect22[xi][yi], rect22[xj][yj]]
           for yi, yj in pairs]
          for xi, xj in pairs]
    elements, crs = check_crs(op, max_generators=4)
    if VERIFY_RULES:
        assert crs.alphabet == "ABCD"
        assert elements == [x + y for x in ("A", "AB", "BA", "B", "") for y in ("C", "CD", "DC", "D", "")]
        assert sorted(crs.rules) == [
            ('AA', 'A'), ('ABA', 'A'),
            ('BAB', 'B'), ('BB', 'B'),
            ('CA', 'AC'), ('CB', 'BC'),
            ('CC', 'C'), ('CDC', 'C'),
            ('DA', 'AD'), ('DB', 'BD'),
            ('DCD', 'D'), ('DD', 'D'),
        ]
        assert crs.essential_counts(5) == [1, 4, 12, 32, 80, 192]
