from mutable_lattice import Vector

from fast_semigroup_homology.cup_products.resolution_with_bar_comparison import (
    ResolutionWithBarComparison
)

def test_rect22():
    op = [[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]]
    res = ResolutionWithBarComparison(op)
    res.extend_to_dimension(10)
    [M0] = res.nodes_by_dimension[0]
    [M1] = res.nodes_by_dimension[1]
    [M2] = res.nodes_by_dimension[2]
    [] = res.nodes_by_dimension[3]

    assert M0.base_node.resolution is res.base_res
    assert M0.base_node.module == [0]
    assert M0.base_node.prev_module is None
    assert M0.base_node.e_images == [Vector([1])]
    assert M0.base_node.children == [M1.base_node]
    assert M0.base_node.child_gen_indexes == [[0]]
    assert M0.e_images_in_bar == [{() : 1}]
    assert M0.prev is None
    assert M0.index_as_child is None

    assert M1.base_node.resolution is res.base_res
    assert M1.base_node.module == [4]
    assert M1.base_node.prev_module == [0]
    assert M1.base_node.e_images == [Vector([1,-1])]
    assert M1.base_node.children == [M2.base_node]
    assert M1.base_node.child_gen_indexes == [[0]]
    assert M1.e_images_in_bar == [{(0,) : 1, (2,): -1}]
    assert M1.prev is M0
    assert M1.index_as_child == 0

    assert M2.base_node.resolution is res.base_res
    assert M2.base_node.module == [0, 0]
    assert M2.base_node.prev_module == [4]
    assert M2.base_node.e_images == [Vector([1,0,0,0,0]), Vector([0,1,0,0,0])]
    assert M2.base_node.children == []
    assert M2.base_node.child_gen_indexes == []
    assert M2.e_images_in_bar == [{(0,0) : 1, (0,2): -1}, {(1,0) : 1, (1,2): -1}]
    assert M2.prev is M1
    assert M2.index_as_child == 0

    assert M0.homology_with_generators_in_bar() == ([0], [{(): 1}])
    assert M1.homology_with_generators_in_bar() == ([], [])
    assert M2.homology_with_generators_in_bar() == ([0], [{(0,0): 1, (0,2): -1, (1,0): -1, (1, 2): 1}])

    assert res.homology_freepart_generators_in_bar(0) == [{(): 1}]
    assert res.homology_freepart_generators_in_bar(1) == []
    assert res.homology_freepart_generators_in_bar(2) == [{(0,0): 1, (0,2): -1, (1,0): -1, (1, 2): 1}]
    assert res.homology_freepart_generators_in_bar(3) == []
    assert res.homology_freepart_generators_in_bar(4) == []
    assert res.homology_freepart_generators_in_bar(5) == []
    assert res.homology_freepart_generators_in_bar(6) == []

def test_C2():
    op = [[0, 1], [1, 0]]
    res = ResolutionWithBarComparison(op)
    res.extend_to_dimension(6)
    [[M0], [M1], [M2], [M3], [M4], [M5], [M6]] = res.nodes_by_dimension
    assert M0.e_images_in_bar == [{(): 1}]
    assert M1.e_images_in_bar == [{(1,): -1}]
    assert M2.e_images_in_bar == [{(1,1): -1}]
    assert M3.e_images_in_bar == [{(1,1,1): 1}]
    assert M4.e_images_in_bar == [{(1,1,1,1): 1}]
    assert M5.e_images_in_bar == [{(1,1,1,1,1): -1}]
    assert M6.e_images_in_bar == [{(1,1,1,1,1,1): -1}]
    assert M0.homology_with_generators_in_bar() == ([0], [{(): 1}])
    assert M1.homology_with_generators_in_bar() == ([2], [{(1,): -1}])
    assert M2.homology_with_generators_in_bar() == ([], [])
    assert M3.homology_with_generators_in_bar() == ([2], [{(1,1,1): 1}])
    assert M4.homology_with_generators_in_bar() == ([], [])
    assert M5.homology_with_generators_in_bar() == ([2], [{(1,1,1,1,1): -1}])
    assert M6.homology_with_generators_in_bar() == ([], [])
    assert res.homology_freepart_generators_in_bar(0) == [{(): 1}]
    assert res.homology_freepart_generators_in_bar(1) == []
    assert res.homology_freepart_generators_in_bar(2) == []
    assert res.homology_freepart_generators_in_bar(3) == []
    assert res.homology_freepart_generators_in_bar(4) == []
    assert res.homology_freepart_generators_in_bar(5) == []
    assert res.homology_freepart_generators_in_bar(7) == []

