from fast_semigroup_homology.cup_products.bar_complex import (
    BarComplex,
)

def test_bar_C2():
    C2 = [[0,1],[1,0]]
    bar = BarComplex(C2)
    assert bar.cup_product_matrix(0, 2) == ([0], [2], [2], [[[1]]])
    assert bar.cup_product_matrix(2, 2) == ([2], [2], [2], [[[1]]])
    assert bar.cup_product_matrix(2, 4) == ([2], [2], [2], [[[1]]])
    assert bar.cup_product_matrix(4, 4) == ([2], [2], [2], [[[1]]])

def test_bar_C3():
    C3 = [[0,1,2],[1,2,0],[2,0,1]]
    bar = BarComplex(C3)
    assert bar.cup_product_matrix(0, 2) == ([0], [3], [3], [[[1]]])
    assert bar.cup_product_matrix(2, 2) == ([3], [3], [3], [[[1]]])
    assert bar.cup_product_matrix(2, 4) == ([3], [3], [3], [[[1]]])
    assert bar.cup_product_matrix(4, 4) == ([3], [3], [3], [[[1]]])

def test_bar_rect22():
    rect22 = [[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]]
    bar = BarComplex(rect22)
    invariants, bar_gens = bar.cohomology_with_generators(0)
    assert invariants == [0]
    assert bar_gens == [{(): 1}]
    invariants, bar_gens = bar.cohomology_with_generators(1)
    assert invariants == []
    assert bar_gens == []
    invariants, bar_gens = bar.cohomology_with_generators(2)
    assert invariants == [0]
    assert bar_gens == [{(0, 0): 1, (0, 1): 1, (2, 0): 1, (2, 1): 1}]

def test_C2():
    res = BarComplex([[0,1],[1,0]])
    assert res.cup_product_exponents(0, 0) == ([0], [0], [[0]])
    assert res.cup_product_exponents(0, 1) == ([0], [], [[]])
    assert res.cup_product_exponents(1, 0) == ([], [0], [])
    assert res.cup_product_exponents(0, 2) == ([0], [2], [[2]])
    assert res.cup_product_exponents(1, 1) == ([], [], [])
    assert res.cup_product_exponents(2, 0) == ([2], [0], [[2]])
    assert res.cup_product_exponents(0, 3) == ([0], [], [[]])
    assert res.cup_product_exponents(1, 2) == ([], [2], [])
    assert res.cup_product_exponents(2, 1) == ([2], [], [[]])
    assert res.cup_product_exponents(3, 0) == ([], [0], [])
    assert res.cup_product_exponents(0, 4) == ([0], [2], [[2]])
    assert res.cup_product_exponents(1, 3) == ([], [], [])
    assert res.cup_product_exponents(2, 2) == ([2], [2], [[2]])
    assert res.cup_product_exponents(3, 1) == ([], [], [])
    assert res.cup_product_exponents(4, 0) == ([2], [0], [[2]])

def test_C3():
    res = BarComplex([[0,1,2],[1,2,0],[2,0,1]])
    assert res.cup_product_exponents(0, 0) == ([0], [0], [[0]])
    assert res.cup_product_exponents(0, 1) == ([0], [], [[]])
    assert res.cup_product_exponents(1, 0) == ([], [0], [])
    assert res.cup_product_exponents(0, 2) == ([0], [3], [[3]])
    assert res.cup_product_exponents(1, 1) == ([], [], [])
    assert res.cup_product_exponents(2, 0) == ([3], [0], [[3]])
    assert res.cup_product_exponents(0, 3) == ([0], [], [[]])
    assert res.cup_product_exponents(1, 2) == ([], [3], [])
    assert res.cup_product_exponents(2, 1) == ([3], [], [[]])
    assert res.cup_product_exponents(3, 0) == ([], [0], [])
    assert res.cup_product_exponents(0, 4) == ([0], [3], [[3]])
    assert res.cup_product_exponents(1, 3) == ([], [], [])
    assert res.cup_product_exponents(2, 2) == ([3], [3], [[3]])
    assert res.cup_product_exponents(3, 1) == ([], [], [])
    assert res.cup_product_exponents(4, 0) == ([3], [0], [[3]])
