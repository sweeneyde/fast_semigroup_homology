from fast_semigroup_homology.cup_products import (
    BarComplex,
    cup_product_freepart_matrix,
    cup_product_matrix,
)

def test_fast_C2():
    C2 = [[0,1],[1,0]]
    assert cup_product_matrix(C2, 0, 0) == ([0], [0], [0], [[[1]]])
    assert cup_product_matrix(C2, 0, 2) == ([0], [2], [2], [[[1]]])
    assert cup_product_matrix(C2, 2, 2) == ([2], [2], [2], [[[1]]])
    for a in [2, 4, 6, 8]:
        assert cup_product_matrix(C2, 0, a) == ([0], [2], [2], [[[1]]])
        assert cup_product_matrix(C2, a, 0) == ([2], [0], [2], [[[1]]])
        for b in [2, 4, 6, 8]:
            assert cup_product_matrix(C2, a, b) == ([2], [2], [2], [[[1]]])

def test_fast_C3():
    C2 = [[0,1,2],[1,2,0],[2,0,1]]
    assert cup_product_matrix(C2, 0, 0) == ([0], [0], [0], [[[1]]])
    assert cup_product_matrix(C2, 0, 2) == ([0], [3], [3], [[[1]]])
    assert cup_product_matrix(C2, 2, 2) == ([3], [3], [3], [[[1]]])
    for a in [2, 4, 6, 8, 10]:
        assert cup_product_matrix(C2, 0, a) == ([0], [3], [3], [[[1]]])
        assert cup_product_matrix(C2, a, 0) == ([3], [0], [3], [[[1]]])
        for b in [2, 4, 6, 8, 10]:
            assert cup_product_matrix(C2, a, b) == ([3], [3], [3], [[[1]]])


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
    for a in [2, 4, 6, 8, 10]:
        assert cup_product_matrix(C3, 0, a) == ([0], [3], [3], [[[1]]])
        assert cup_product_matrix(C3, a, 0) == ([3], [0], [3], [[[1]]])
        for b in [2, 4, 6, 8, 10]:
            assert cup_product_matrix(C3, a, b) == ([3], [3], [3], [[[1]]])

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

def test_rect22_rect22_slow():
    rect22 = [[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]]
    elements = [(i, j) for i in range(5) for j in range(5)]
    element_to_index = {ij: index for index, ij in enumerate(elements)}
    op = [[element_to_index[rect22[xi][yi], rect22[xj][yj]]
           for yi, yj in elements]
          for xi, xj in elements]
    mat = cup_product_freepart_matrix(op, 2, 2)
    [[[a],[b]],[[c],[d]]] = mat
    # For invertible integer change of basis matrix P,
    # we have abs(det(P)) == 1, so upon changing basis, we have
    #   det(P * A * P^T) == det(P) * det(A) * det(P) == det(A).
    # Changing basis didn't change the determinant.
    # In the standard basis of generators for the homology
    # of S^2 x S^2, we would have the matrix [[0,1],[1,0]],
    # so check that we at least have the right determinant:
    assert a*d - b*c == -1, (a,b,c,d)

def test_rect22_rect22_fast():
    rect22 = [[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]]
    elements = [(i, j) for i in range(5) for j in range(5)]
    element_to_index = {ij: index for index, ij in enumerate(elements)}
    op = [[element_to_index[rect22[xi][yi], rect22[xj][yj]]
           for yi, yj in elements]
          for xi, xj in elements]
    inva, invb, invc, mat = cup_product_matrix(op, 2, 2, max_generators=4)
    assert inva == invb == [0, 0]
    assert invc == [0]
    [[[a],[b]],[[c],[d]]] = mat
    assert a*d - b*c == -1, (a,b,c,d)

# def test_Z0ZZZ_slow():
#     op = [[0,1,0,1,0,0],
#           [0,1,0,1,0,1],
#           [2,3,2,3,2,2],
#           [2,3,2,3,2,3],
#           [0,1,0,1,0,4],
#           [0,1,2,3,4,5]]
#     [[[c02]]] = cup_product_freepart_matrix(op, 0, 2)
#     assert abs(c02) == 1
#     [[[c03]]] = cup_product_freepart_matrix(op, 0, 3)
#     assert abs(c03) == 1
#     for a in [2,3,4]:
#         for b in [2,3,4]:
#             assert cup_product_freepart_matrix(op, a, b) == [[[0]]]

def test_Z0ZZZ_fast():
    op = [[0,1,0,1,0,0],
          [0,1,0,1,0,1],
          [2,3,2,3,2,2],
          [2,3,2,3,2,3],
          [0,1,0,1,0,4],
          [0,1,2,3,4,5]]
    for a in [2,3,4]:
        for b in [2,3,4]:
            assert cup_product_matrix(op, a, b) == ([0], [0], [0], [[[0]]])

def test_rect22_adjoin_C2_2sided_action():
    op = [[0,1,0,1,0,1],
          [0,1,0,1,1,0],
          [2,3,2,3,2,3],
          [2,3,2,3,3,2],
          [0,1,2,3,4,5],
          [2,3,0,1,5,4]]
    assert cup_product_matrix(op, 2, 2) == ([0], [0], [2], [[[1]]])
    assert cup_product_matrix(op, 2, 4) == ([0], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 2, 6) == ([0], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 4, 4) == ([2], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 2, 8) == ([0], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 4, 6) == ([2], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 2, 10) == ([0], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 4, 8) == ([2], [2], [2], [[[1]]])
    assert cup_product_matrix(op, 6, 6) == ([2], [2], [2], [[[1]]])

def test_rect22_adjoin_C2_trivially():
    op = [[0,1,0,1,0,0],
          [0,1,0,1,1,1],
          [2,3,2,3,2,2],
          [2,3,2,3,3,3],
          [0,1,2,3,4,5],
          [0,1,2,3,5,4]]
    assert cup_product_matrix(op, 2, 2) == ([0], [0], [2], [[[0]]])
    assert cup_product_matrix(op, 2, 4) == ([0], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 2, 6) == ([0], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 4, 4) == ([2], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 2, 8) == ([0], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 4, 6) == ([2], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 2, 10) == ([0], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 4, 8) == ([2], [2], [2], [[[0]]])
    # assert cup_product_matrix(op, 6, 6) == ([2], [2], [2], [[[0]]])



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

# def test_SusC2():
#     res = BarComplex([[0,1,0,1,0,1],[0,1,0,1,1,0],[2,3,2,3,2,3],[2,3,2,3,3,2],[0,1,2,3,4,5],[0,1,2,3,5,4]])
#     assert res.cup_product_exponents(3, 3) == ([2], [2], [[1]])