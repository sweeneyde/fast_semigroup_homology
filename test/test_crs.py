from fast_semigroup_homology.cup_products.rewriting.crs import CRS
from itertools import product

def check_crs_consistency(crs):
    essentials = set()
    for cells in crs._essentials:
        essentials.update(cells)
    [root] = crs._essentials[0]
    stack = [root]
    cells = set()
    while stack:
        c = stack.pop()
        if c is None or c in cells:
            continue
        cells.add(c)
        stack.extend(c._extensions.values())
    assert essentials <= cells
    tp_root = type(root)
    for cell in cells:
        for ext in cell._extensions.values():
            assert ext.all_but_last is cell
        if cell is root:
            continue
        abl = cell.all_but_last
        assert cell is abl._extensions[cell.last]
        if type(cell) is tp_root:
            assert cell in essentials
            assert type(abl) is tp_root
            assert len(cell.faces) == len(abl.faces) + 1
            assert cell.faces[-1] is abl
        else:
            assert cell not in essentials
            assert cell.paired_coface_faces.count(cell) == 1
            if type(abl) is tp_root:
                assert len(cell.paired_coface_faces) == len(abl.faces) + 2
            else:
                assert len(cell.paired_coface_faces) == len(abl.paired_coface_faces) + 1

def test_trivial():
    crs = CRS("", {})
    assert list(crs.elements()) == [""]
    assert crs.essential_tuples(4) == [[()], [], [], [], []]
    assert crs.essential_counts(4) == [1, 0, 0, 0, 0]
    assert crs.op_table([""]) == [[0]]
    assert crs.homology_list(5) == [{0: 1}, {}, {}, {}, {}, {}]
    assert crs.cohomology_list(5) == [{0: 1}, {}, {}, {}, {}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    check_crs_consistency(crs)

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
    check_crs_consistency(crs)
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
    check_crs_consistency(crs)
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
    check_crs_consistency(crs)

    [[root], [x], [x_xx], [x_xx_x], [x_xx_x_xx]] = crs._essentials[:5]
    xx = root._extensions["xx"]
    xx_x = xx._extensions["x"]
    xx_x_xx = xx_x._extensions["xx"]
    x_xx_xx = x_xx._extensions["xx"]
    xx_x_x = xx_x._extensions["x"]
    assert root.faces == (None,)
    assert x.faces == (root, root)
    assert x_xx.faces == (xx, None, x)
    assert x_xx_x.faces == (xx_x, None, None, x_xx)
    assert x_xx_x_xx.faces == (xx_x_xx, None, None, None, x_xx_x)
    assert xx.paired_coface_faces == (x, xx, x)
    assert xx_x.paired_coface_faces == (None, xx_x, x_xx, None)
    assert xx_x_xx.paired_coface_faces == (None, xx_x_xx, x_xx_xx, None, None)
    assert x_xx_xx.paired_coface_faces == (xx_x_x, None, None, x_xx_xx, x_xx_x)
    assert xx_x_x.paired_coface_faces == (None, xx_x_x, x_xx_x, None, None)

    assert crs.essential_representation(("xx",)) == {("x",): 2}
    assert crs.essential_representation(("xx", "x")) == {("x", "xx"): 1}
    assert crs.essential_representation(("xx", "x", "x")) == {("x", "xx", "x"): 1}
    assert crs.essential_representation(("x", "xx", "xx")) == {("x", "xx", "x"): 2}
    assert crs.essential_representation(("xx", "x", "xx")) == {("x", "xx", "x"): 2}
    for e in [("x",), ("x","xx"), ("x","xx","x"), ("x","xx","x","xx")]:
        assert crs.essential_representation(e) == {e: 1}
    assert crs.boundary_dicts(0) == [{}]
    assert crs.boundary_dicts(1) == [{}]
    assert crs.boundary_dicts(2) == [{0: 3}]
    assert crs.boundary_dicts(3) == [{}]

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
    check_crs_consistency(crs)

def test_right_zero():
    crs = CRS("xy", {"xx":"x", "xy":"y", "yx":"x", "yy":"y"})
    assert list(crs.elements()) == ["", "x", "y"]
    assert crs.essential_tuples(4) == [list(product("xy", repeat=n)) for n in range(4+1)]
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert crs.op_table(["", "x", "y"]) == [[0, 1, 2], [1, 1, 2], [2, 1, 2]]
    assert crs.homology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.cohomology_list(3) == [{0: 1}, {}, {}, {}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    check_crs_consistency(crs)

def test_CRS_rect22():
    crs = CRS("xy", {"xx":"x", "yy":"y", "xyx":"x", "yxy":"y"})
    assert crs.essential_counts(4) == [1, 2, 4, 8, 16]
    assert crs.essential_tuples(3) == [
        [()],
        [("x",), ("y",)],
        [("x", "x"), ("x", "yx"), ("y", "y"), ("y", "xy")],
        [("x","x","x"), ("x","x","yx"), ("x","yx","y"), ("x","yx","x"),
         ("y","y","y"), ("y","y","xy"), ("y","xy","x"), ("y","xy","y")]
    ]
    check_crs_consistency(crs)

    [[root], [x,y], [x_x,x_yx,y_y,y_xy],
     [x_x_x,x_x_yx,x_yx_y,x_yx_x,y_y_y,y_y_xy,y_xy_x,y_xy_y]] = crs._essentials[:4]
    assert root.faces == (None,)
    assert x.faces == (root, root)
    assert y.faces == (root, root)
    assert x_x.faces == (x, x, x)
    xy = root._extensions["xy"]
    yx = root._extensions["yx"]
    assert xy.paired_coface_faces == (y, xy, x)
    assert yx.paired_coface_faces == (x, yx, y)
    assert x_yx.faces == (yx, x, x)
    assert y_y.faces == (y, y, y)
    assert y_xy.faces == (xy, y, y)
    assert x_x_x.faces == (x_x, x_x, x_x, x_x)
    assert x_x_yx.faces == (x_yx, x_yx, x_x, x_x)
    xy_x = xy._extensions["x"]
    xy_y = xy._extensions["y"]
    yx_y = yx._extensions["y"]
    yx_x = yx._extensions["x"]
    assert xy_x.paired_coface_faces == (None, xy_x, x_yx, None)
    assert xy_y.paired_coface_faces == (y_y, xy_y, None, None)
    assert yx_y.paired_coface_faces == (None, yx_y, y_xy, None)
    assert yx_x.paired_coface_faces == (x_x, yx_x, None, None)
    assert x_yx_y.faces == (yx_y, None, None, x_yx)
    assert x_yx_x.faces == (yx_x, x_x, x_yx, x_yx)
    assert y_y_y.faces == (y_y, y_y, y_y, y_y)
    assert y_y_xy.faces == (y_xy, y_xy, y_y, y_y)
    assert y_xy_x.faces == (xy_x, None, None, y_xy)
    assert y_xy_y.faces == (xy_y, y_y, y_xy, y_xy)

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
    check_crs_consistency(crs)
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
    check_crs_consistency(crs)
    assert crs.essential_counts(4) == [1, 3, 7, 15, 31]
    assert crs.homology_list(3) == [{0: 1}, {}, {2: 1}, {}]
    assert crs.cohomology_list(3) == [{0: 1}, {}, {}, {2: 1}]
    assert crs.cohomology_with_generators(0) == ([0], [{(): 1}])
    assert crs.cohomology_with_generators(3) == ([2], [{("g", "g", "y"): 1}])

def test_infinite_example():
    crs = CRS('xy', {'xxx':'', 'yyy':'', 'yxyx':'xxyy', "yyxx":"xyxy"})
    crs.compute_essentials(5)
    assert crs.essential_tuples(3) == [
        [()],
        [("x",), ("y",)],
        [("x", "xx"), ("y", "yy"), ("y", "xyx"), ("y", "yxx")],
        [("x", "xx", "x"), ("y", "yy", "y"),
         ("y", "yy", "xx"), ("y", "yy", "xyx"),
         ("y", "xyx", "yx"), ("y", "xyx", "xx"), ("y", "yxx", "x")]
    ]
    check_crs_consistency(crs)
