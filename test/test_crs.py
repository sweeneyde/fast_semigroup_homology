from fast_semigroup_homology.cup_products.rewriting.crs import CRS
from itertools import product

def test_trivial():
    crs = CRS("", {})
    assert list(crs.elements()) == [""]
    assert crs.essential_tuples(4) == [[()], [], [], [], []]
    assert crs.essential_counts(4) == [1, 0, 0, 0, 0]
    assert crs.op_table([""]) == [[0]]
    assert crs.homology_list(5) == [{0: 1}, {}, {}, {}, {}, {}]
    assert crs.cohomology_list(5) == [{0: 1}, {}, {}, {}, {}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])

def test_onegen():
    crs = CRS("x", {"xx":"x"})
    assert list(crs.elements()) == ["", "x"]
    assert crs.essential_tuples(4) == [
        [()],
        [("x",)],
        [("x","x")],
        [("x","x","x")],
        [("x","x","x","x")],
    ]
    assert crs.essential_counts(4) == [1, 1, 1, 1, 1]
    assert crs.op_table(["", "x"]) == [[0, 1], [1, 1]]
    assert crs.op_table(["x", ""]) == [[0, 0], [0, 1]]
    assert crs.homology_list(5) == [{0: 1}, {}, {}, {}, {}, {}]
    assert crs.cohomology_list(5) == [{0: 1}, {}, {}, {}, {}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])

def test_C2():
    crs = CRS("x", {"xx":""})
    assert list(crs.elements()) == ["", "x"]
    assert crs.essential_tuples(4) == [
        [()],
        [("x",)],
        [("x","x")],
        [("x","x","x")],
        [("x","x","x","x")],
    ]
    assert crs.essential_counts(4) == [1, 1, 1, 1, 1]
    assert crs.boundary_nonzero_invariants(0) == []
    assert crs.boundary_nonzero_invariants(1) == []
    assert crs.boundary_nonzero_invariants(2) == [2]
    assert crs.boundary_nonzero_invariants(3) == []
    assert crs.boundary_nonzero_invariants(4) == [2]
    assert crs.op_table(["", "x"]) == [[0, 1], [1, 0]]
    assert crs.homology_list(5) == [{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}]
    assert crs.cohomology_list(5) == [{0: 1}, {}, {2: 1}, {}, {2: 1}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    assert crs.cohomology_with_generators(2) == ([2], [{("x", "x"): 1}])
    assert crs.cohomology_with_generators(4) == ([2], [{("x", "x", "x", "x"): 1}])

def test_C3():
    crs = CRS("x", {"xxx":""})
    assert list(crs.elements()) == ["", "x", "xx"]
    assert crs.essential_tuples(4) == [
        [()],
        [("x",)],
        [("x","xx")],
        [("x","xx","x")],
        [("x","xx","x","xx")],
    ]
    assert crs.essential_counts(4) == [1, 1, 1, 1, 1]
    assert crs.op_table(["", "x", "xx"]) == [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
    [[root], [x], [x_xx], [x_xx_x], [x_xx_x_xx]] = crs._essentials[:5]
    assert root._extensions == {"x": x}
    assert x._extensions == {"xx": x_xx}
    assert x_xx._extensions == {"x": x_xx_x}
    assert x_xx_x._extensions == {"xx": x_xx_x_xx}
    assert x.last == x_xx_x.last == "x"
    assert x_xx.last is x_xx_x_xx.last == "xx"
    assert x.get_faces(crs) == (root, root)
    xx = root.get_extension("xx", crs)
    assert x_xx.get_faces(crs) == (xx, None, x)
    assert xx._essrep == {x: 2}
    assert xx._inessrep == {}
    assert crs.homology_list(5) == [{0: 1}, {3: 1}, {}, {3: 1}, {}, {3: 1}]
    assert crs.cohomology_list(5) == [{0: 1}, {}, {3: 1}, {}, {3: 1}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    assert crs.cohomology_with_generators(2) == ([3], [{("x", "xx"): 1}])
    assert crs.cohomology_with_generators(4) == ([3], [{("x", "xx", "x", "xx"): 1}])

def test_left_zero():
    crs = CRS("xy", {"xx":"x", "xy":"x", "yx":"y", "yy":"y"})
    assert list(crs.elements()) == ["", "x", "y"]
    assert crs.essential_tuples(4) == [list(product("xy", repeat=n)) for n in range(4+1)]
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert crs.op_table(["", "x", "y"]) == [[0, 1, 2], [1, 1, 1], [2, 2, 2]]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.cohomology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])

def test_right_zero():
    crs = CRS("xy", {"xx":"x", "xy":"y", "yx":"x", "yy":"y"})
    assert list(crs.elements()) == ["", "x", "y"]
    assert crs.essential_tuples(4) == [list(product("xy", repeat=n)) for n in range(4+1)]
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert crs.op_table(["", "x", "y"]) == [[0, 1, 2], [1, 1, 2], [2, 1, 2]]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.cohomology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])

def test_CRS_rect22():
    crs = CRS("xy", {"xx":"x", "yy":"y", "xyx":"x", "yxy":"y"})
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert list(crs.elements()) == ["", "x", "y", "xy", "yx"]
    assert crs.homology_list(3) == [{0: 1}, {}, {0: 1}, {}]
    assert crs.op_table(["", "x", "y", "xy", "yx"]) == (
        [[0, 1,2,3,4], [1, 1,3,3,1], [2, 4,2,2,4], [3, 1,3,3,1], [4, 4,2,2,4]]
    )
    assert crs.op_table(["x", "xy", "yx", "y", ""]) == (
        [[0,1,0,1, 0], [0,1,0,1, 1], [2,3,2,3, 2], [2,3,2,3, 3], [0,1,2,3, 4]]
    )
    assert crs.homology_list(3) == [{0: 1}, {}, {0: 1}, {}]
    assert crs.cohomology_list(3) == [{0: 1}, {}, {0: 1}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    assert crs.cohomology_with_generators(2) == ([0], [{("x", "x"): 1}])

def test_rect22_with_C2_acting_on_both_sides():
    crs = CRS("gx", {"xx":"x", "gg":"", "xgx":"x"})
    assert list(crs.elements()) == ["", "g", "x", "gx", "xg", "gxg"]
    assert crs.essential_tuples(4) == [
        [()],
        [("g",), ("x",)],
        [("g","g"), ("x","x"), ("x", "gx")],
        [("g","g","g"), ("x","x","x"), ("x","x","gx"), ("x","gx","x"), ("x","gx","gx")],
        [("g","g","g","g"),
         ("x","x","x","x"), ("x","x","x","gx"), ("x","x","gx","x"), ("x","x","gx","gx"),
         ("x","gx","x","x"), ("x","gx","x","gx"), ("x","gx","gx","x"), ("x","gx","gx","gx")],
    ]
    assert crs.essential_counts(4) == [1, 2, 3, 5, 9]
    assert crs.op_table(["x", "xg", "gx", "gxg", "", "g"]) == (
        [[0,1,0,1,0,1],
         [0,1,0,1,1,0],
         [2,3,2,3,2,3],
         [2,3,2,3,3,2],
         [0,1,2,3,4,5],
         [2,3,0,1,5,4]]
    )
    assert crs.homology_list(5) == [{0: 1}, {}, {0: 1}, {2: 1}, {}, {2: 1}]
    assert crs.cohomology_list(5) == [{0: 1}, {}, {0: 1}, {}, {2: 1}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    assert crs.cohomology_with_generators(2) == ([0], [{("g", "g"): 1}])
    assert crs.cohomology_with_generators(4) == ([2], [{("g", "g", "g", "g"): 1}])

def test_rect22_with_C2_acting_on_one_side():
    crs = CRS("gxy", {"gg":"", "gx":"x", "gy":"y",
                      "xx":"x", "xy":"x", "yx":"y", "yy":"y"})
    assert list(crs.elements()) == ["", "g", "x", "y", "xg", "yg"]
    assert crs.essential_tuples(3) == [
        [()],
        [("g",), ("x",), ("y",)],
        [("g","g"), ("g","x"), ("g","y"), ("x","x"), ("x","y"), ("y","x"), ("y","y")],
        [("g","g","g"), ("g","g","x"), ("g","g","y"),
         ("g","x","x"), ("g","x","y"), ("g","y","x"), ("g","y","y"),
         ("x","x","x"), ("x","x","y"), ("x","y","x"), ("x","y","y"),
         ("y","x","x"), ("y","x","y"), ("y","y","x"), ("y","y","y")]
    ]
    assert crs.essential_counts(4) == [1, 3, 7, 15, 31]
    assert crs.homology_list(3) == [{0: 1}, {}, {2: 1}, {}]
    assert crs.cohomology_list(3) == [{0: 1}, {}, {}, {2: 1}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    assert crs.cohomology_with_generators(3) == ([2], [{("g", "g", "y"): 1}])
