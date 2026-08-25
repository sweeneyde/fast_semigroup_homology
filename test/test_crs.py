from fast_semigroup_homology.cup_products.rewriting.crs import CRS
from itertools import product

def test_trivial():
    crs = CRS("", [])
    assert list(crs.elements()) == [""]
    crs.compute_essentials(4)
    assert crs.essentials == [((),), (), (), (), ()]
    assert crs.essential_counts(4) == [1, 0, 0, 0, 0]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.op_table_with_map() == (
        [""],
        {"": 0},
        [[0]]
    )

def test_onegen():
    crs = CRS("x", [("xx", "x")])
    assert list(crs.elements()) == ["", "x"]
    crs.compute_essentials(4)
    assert crs.essentials == [
        ((),),
        (("x",),),
        (("x","x"),),
        (("x","x","x"),),
        (("x","x","x","x"),),
    ]
    assert crs.essential_counts(4) == [1, 1, 1, 1, 1]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.op_table_with_map() == (
        ["", "x"],
        {"": 0, "x": 1},
        [[0, 1], [1, 1]],
    )
    assert crs.op_table_with_map(["x", ""]) == (
        ["x", ""],
        {"x": 0, "": 1},
        [[0, 0], [0, 1]],
    )

def test_C2():
    crs = CRS("x", [("xx", "")])
    assert list(crs.elements()) == ["", "x"]
    crs.compute_essentials(4)
    assert crs.essentials == [
        ((),),
        (("x",),),
        (("x","x"),),
        (("x","x","x"),),
        (("x","x","x","x"),),
    ]
    assert crs.essential_counts(4) == [1, 1, 1, 1, 1]
    assert crs.boundary_nonzero_invariants(0) == []
    assert crs.boundary_nonzero_invariants(1) == []
    assert crs.boundary_nonzero_invariants(2) == [2]
    assert crs.boundary_nonzero_invariants(3) == []
    assert crs.boundary_nonzero_invariants(4) == [2]
    assert crs.homology_list(3) == [{0: 1}, {2: 1}, {}, {2: 1}]
    assert crs.op_table_with_map() == (
        ["", "x"],
        {"": 0, "x": 1},
        [[0, 1], [1, 0]],
    )

def test_C3():
    crs = CRS("x", [("xxx", "")])
    assert list(crs.elements()) == ["", "x", "xx"]
    crs.compute_essentials(4)
    assert crs.essentials == [
        ((),),
        (("x",),),
        (("x","xx"),),
        (("x","xx","x"),),
        (("x","xx","x","xx"),),
    ]
    assert crs.essential_counts(4) == [1, 1, 1, 1, 1]
    assert crs.homology_list(3) == [{0: 1}, {3: 1}, {}, {3: 1}]
    assert crs.op_table_with_map() == (
        ["", "x", "xx"],
        {"": 0, "x": 1, "xx": 2},
        [[0, 1, 2], [1, 2, 0], [2, 0, 1]],
    )


def test_left_zero():
    crs = CRS("xy", [("xx", "x"), ("xy", "x"), ("yx", "y"), ("yy", "y")])
    assert list(crs.elements()) == ["", "x", "y"]
    crs.compute_essentials(4)
    assert crs.essentials == [tuple(product("xy", repeat=n)) for n in range(4+1)]
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.op_table_with_map() == (
        ["", "x", "y"],
        {"": 0, "x": 1, "y": 2},
        [[0, 1, 2], [1, 1, 1], [2, 2, 2]],
    )

def test_right_zero():
    crs = CRS("xy", [("xx", "x"), ("xy", "y"), ("yx", "x"), ("yy", "y")])
    assert list(crs.elements()) == ["", "x", "y"]
    crs.compute_essentials(4)
    assert crs.essentials == [tuple(product("xy", repeat=n)) for n in range(4+1)]
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.op_table_with_map() == (
        ["", "x", "y"],
        {"": 0, "x": 1, "y": 2},
        [[0, 1, 2], [1, 1, 2], [2, 1, 2]],
    )

def test_CRS_rect22():
    crs = CRS("xy", [("xx", "x"), ("yy", "y"), ("xyx", "x"), ("yxy", "y")])
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert list(crs.elements()) == ["", "x", "y", "xy", "yx"]
    assert crs.homology_list(3) == [{0: 1}, {}, {0: 1}, {}]
    assert crs.op_table_with_map() == (
        ["", "x", "y", "xy", "yx"],
        {"": 0, "x": 1, "y": 2, "xy": 3, "yx": 4},
        [[0, 1,2,3,4], [1, 1,3,3,1], [2, 4,2,2,4], [3, 1,3,3,1], [4, 4,2,2,4]],
    )
    assert crs.op_table_with_map(["x", "xy", "yx", "y", ""]) == (
        ["x", "xy", "yx", "y", ""],
        {"x": 0, "xy": 1, "yx": 2, "y": 3, "": 4},
        [[0,1,0,1, 0], [0,1,0,1, 1], [2,3,2,3, 2], [2,3,2,3, 3], [0,1,2,3, 4]],
    )

def test_rect22_with_C2_acting_on_both_sides():
    crs = CRS("gx", [("xx", "x"), ("gg", ""), ("xgx", "x")])
    assert list(crs.elements()) == ["", "g", "x", "gx", "xg", "gxg"]
    crs.compute_essentials(4)
    assert crs.essentials == [
        ((),),
        (("g",), ("x",)),
        (("g","g"), ("x","x"), ("x", "gx")),
        (("g","g","g"), ("x","x","x"), ("x","x","gx"), ("x","gx","x"), ("x","gx","gx")),
        (("g","g","g","g"),
         ("x","x","x","x"), ("x","x","x","gx"), ("x","x","gx","x"), ("x","x","gx","gx"),
         ("x","gx","x","x"), ("x","gx","x","gx"), ("x","gx","gx","x"), ("x","gx","gx","gx")),
    ]
    assert crs.essential_counts(4) == [1, 2, 3, 5, 9]
    assert crs.homology_list(3) == [{0: 1}, {}, {0: 1}, {2: 1}]
    assert crs.op_table_with_map(["x", "xg", "gx", "gxg", "", "g"]) == (
        ["x", "xg", "gx", "gxg", "", "g"],
        {"x": 0, "xg": 1, "gx": 2, "gxg": 3, "": 4, "g": 5},
        [[0,1,0,1,0,1],
         [0,1,0,1,1,0],
         [2,3,2,3,2,3],
         [2,3,2,3,3,2],
         [0,1,2,3,4,5],
         [2,3,0,1,5,4]]
    )

def test_rect22_with_C2_acting_on_one_side():
    crs = CRS("gxy", [("gg", ""), ("gx", "x"), ("gy", "y"),
                      ("xx", "x"), ("xy", "x"), ("yx", "y"), ("yy", "y")])
    assert list(crs.elements()) == ["", "g", "x", "y", "xg", "yg"]
    crs.compute_essentials(3)
    assert crs.essentials == [
        ((),),
        (("g",), ("x",), ("y",)),
        (("g","g"), ("g","x"), ("g","y"), ("x","x"), ("x","y"), ("y","x"), ("y","y")),
        (("g","g","g"), ("g","g","x"), ("g","g","y"),
         ("g","x","x"), ("g","x","y"), ("g","y","x"), ("g","y","y"),
         ("x","x","x"), ("x","x","y"), ("x","y","x"), ("x","y","y"),
         ("y","x","x"), ("y","x","y"), ("y","y","x"), ("y","y","y"),),
    ]
    assert crs.essential_counts(4) == [1, 3, 7, 15, 31]
    assert crs.homology_list(3) == [{0: 1}, {}, {2: 1}, {}]
