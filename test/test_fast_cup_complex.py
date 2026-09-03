from fast_semigroup_homology.cup_products import FastCupComplex as Cup
from fast_semigroup_homology.cup_products.rewriting.crs import CRS

import pytest

@pytest.mark.parametrize("n", [2,3,4])
def test_cyclic(n):
    c = Cup(CRS("x", [("x"*n, "")]))
    assert c.op == [[(i+j)%n for j in range(n)] for i in range(n)]
    assert c.cohomology_list(10) == [[0],[],[n],[],[n],[],[n],[],[n],[],[n]]
    for a in range(10):
        for b in range(10):
            matrix = c.cup_product_matrix(a, b)
            if a % 2 == 1:
                assert matrix == []
            elif b % 2 == 1:
                assert matrix == [[]]
            else:
                assert matrix == [[[1]]]

def test_V4():
    c = Cup(CRS("xy", [("yx", "xy"), ("xx", ""), ("yy", "")]))
    assert c.op == [[i ^ j for j in range(4)] for i in range(4)]
    assert c.cohomology_list(6) == [[0],[],[2,2],[2],[2,2,2],[2,2],[2,2,2,2]]
    [[aa, ab], [ba, bb]] = c.cup_product_matrix(2, 2)
    assert ab == ba
    # linearly independent over Z/2
    assert ab != [x^y for x, y in zip(aa, bb)]
    assert any(aa) and any(ab) and any(bb)
    assert aa != ab and aa != bb and ab != bb

def test_rect22_rect22():
    # Product of (2x2 rectangular band + 1) with itself
    c = Cup(CRS("abxy", [("aa","a"),("bb","b"),("aba","a"),("bab","b"),
                         ("xx","x"),("yy","y"),("xyx","x"),("yxy","y"),
                         ("xa","ax"),("xb","bx"),("ya","ay"),("yb","by")]))
    assert len(c.op) == 25
    assert c.cohomology_list(10) == [[0],[],[0,0],[],[0],[],[],[],[],[],[]]
    assert c.cup_product_matrix(0, 0) == [[[1]]]
    assert c.cup_product_matrix(0, 1) == [[]]
    assert c.cup_product_matrix(0, 2) == [[[1,0],[0,1]]]
    assert c.cup_product_matrix(0, 3) == [[]]
    assert c.cup_product_matrix(0, 4) == [[[-1]]] # opposite generator is okay
    assert c.cup_product_matrix(0, 5) == [[]]
    for b in range(6):
        assert c.cup_product_matrix(1, b) == []
        assert c.cup_product_matrix(3, b) == []
    assert c.cup_product_matrix(2, 0) == [[[1,0]],[[0,1]]]
    assert c.cup_product_matrix(2, 1) == [[], []]
    [[[x],[y]],[[z],[w]]] = c.cup_product_matrix(2, 2)
    assert x*w - y*z == -1, (x,y,z,w)
    assert c.cup_product_matrix(2, 3) == [[], []]
    assert c.cup_product_matrix(2, 4) == [[[]], [[]]]
    assert c.cup_product_matrix(2, 5) == [[], []]

def test_SusC2():
    # Suspensions have no nontrivial cup producsts
    c = Cup([[0,1,0,1,0,1],
             [0,1,0,1,1,0],
             [2,3,2,3,2,3],
             [2,3,2,3,3,2],
             [0,1,2,3,4,5],
             [0,1,2,3,5,4]])
    assert c.cohomology_list(10) == [[0],[],[],[2],[],[2],[],[2],[],[2],[]]
    for a in range(1, 5):
        for b in range(1, 5):
            matrix = c.cup_product_matrix(a, b)
            assert not any(any(cell) for row in matrix for cell in row)
    assert c.cup_product_matrix(0, 0) == [[[1]]]
    assert c.cup_product_matrix(0, 1) == [[]]
    assert c.cup_product_matrix(0, 2) == [[]]
    assert c.cup_product_matrix(0, 3) == [[[1]]]
    assert c.cup_product_matrix(0, 4) == [[]]
    assert c.cup_product_matrix(0, 5) == [[[1]]]
    assert c.cup_product_matrix(1, 0) == []
    assert c.cup_product_matrix(2, 0) == []
    assert c.cup_product_matrix(3, 0) == [[[1]]]
    assert c.cup_product_matrix(4, 0) == []
    assert c.cup_product_matrix(5, 0) == [[[1]]]

def test_simple_2_2_C2_sandwich0():
    # See Steinberg: https://arxiv.org/pdf/2405.06594
    # This simple semigroup is B(C2) v S2
    c = Cup([[0,1,2,3,0,1,2,3,0],
             [1,0,3,2,1,0,3,2,1],
             [0,1,2,3,0,1,2,3,2],
             [1,0,3,2,1,0,3,2,3],
             [4,5,6,7,4,5,6,7,4],
             [5,4,7,6,5,4,7,6,5],
             [4,5,6,7,4,5,6,7,6],
             [5,4,7,6,5,4,7,6,7],
             [0,1,2,3,4,5,6,7,8]], max_generators=2)
    assert c.cohomology_list(10) == [[0],[],[2,0],[],[2],[],[2],[],[2],[],[2]]
    # The free generator of H^2 could correspond to S2...
    assert c.cup_product_matrix(2, 2) == [[[1],[0]],[[0],[0]]]
    assert c.cup_product_matrix(2, 4) == [[[1]],[[0]]]
    # ...or could be a sum of that generator and the torsion element
    # assert c.cup_product_matrix(2, 2) == [[[1],[1]],[[1],[1]]]
    # assert c.cup_product_matrix(2, 4) == [[[1]],[[1]]]

def test_simple_2_2_C2_sandwich1():
    # See Steinberg: https://arxiv.org/pdf/2405.06594
    # This simple semigroup is B(C2) with a 2-disk glued in along the generator of pi_1.
    c = Cup([[0,1,2,3,0,1,2,3,0],
             [1,0,3,2,1,0,3,2,1],
             [0,1,2,3,1,0,3,2,2],
             [1,0,3,2,0,1,2,3,3],
             [4,5,6,7,4,5,6,7,4],
             [5,4,7,6,5,4,7,6,5],
             [4,5,6,7,5,4,7,6,6],
             [5,4,7,6,4,5,6,7,7],
             [0,1,2,3,4,5,6,7,8]], max_generators=2)
    assert c.cohomology_list(10) == [[0],[],[0],[],[2],[],[2],[],[2],[],[2]]
    assert c.cup_product_matrix(2, 2) == [[[1]]]
    assert c.cup_product_matrix(2, 4) == [[[1]]]
    assert c.cup_product_matrix(2, 6) == [[[1]]]
    assert c.cup_product_matrix(4, 4) == [[[1]]]
