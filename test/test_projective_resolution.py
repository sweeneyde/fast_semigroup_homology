from fast_monoid_homology.projective_resolution import ProjectiveResolution
from mutable_lattice import Vector

def test_trivial():
    res = ProjectiveResolution([[0]])
    assert res.homology_list(10) == [{0: 1}] + [{}] * 10
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    [M1] = M0.children
    assert M1.module == []
    assert M1.prev_module == [0]
    assert M1.e_images == []
    assert M1.children == []

def test_C2():
    res = ProjectiveResolution([[0,1],[1,0]])
    assert res.homology_list(10) == [{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}]
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    [M1] = M0.children
    assert M1.module == [0]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1, -1])]
    [M2] = M1.children
    assert M2.module == [0]
    assert M2.prev_module == [0]
    assert M2.e_images == [Vector([1, 1])]
    [M3] = M2.children
    assert M3 is M1

def test_C3():
    res = ProjectiveResolution([[0,1,2],[1,2,0],[2,0,1]])
    assert res.homology_list(10) == [{0: 1}] + [{3: 1}, {}] * 5
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    [M1] = M0.children
    assert M1.module == [0]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1,0,-1])]
    [M2] = M1.children
    assert M2.module == [0]
    assert M2.prev_module == [0]
    assert M2.e_images == [Vector([1, 1, 1])]
    [M3] = M2.children
    assert M3 is M1

def test_rect22():
    res = ProjectiveResolution([[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]])
    assert res.homology_list(10) == [{0: 1}, {}, {0: 1}, {}, {}, {}, {}, {}, {}, {}, {}]
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    [M1] = M0.children
    assert M1.module == [4]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1, -1])]
    [M2] = M1.children
    assert M2.module == [0, 0]
    assert M2.prev_module == [4]
    assert M2.e_images == [Vector([1,0,0,0,0]),
                           Vector([0,1,0,0,0])]
    [M3a, M3b] = M2.children
    assert M3a is M3b
    assert M3a.module == []
    assert M3a.prev_module == [0]
    assert M3a.e_images == []
    assert M3a.children == []

def test_infinitely_many_Zs():
    res = ProjectiveResolution([[0,1,0,1,0,0],[0,1,0,1,0,1],[2,3,2,3,2,2],[2,3,2,3,2,3],[0,1,0,1,0,4],[0,1,2,3,4,5]])
    assert res.homology_list(10) == [{0: 1}, {}] + [{0: 1}] * 9
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    [M1] = M0.children
    assert M1.module == [5]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1,-1])]
    [M2] = M1.children
    assert M2.module == [0, 5]
    assert M2.prev_module == [5]
    assert M2.e_images == [Vector([0,1,0,0,0,0]),
                           Vector([0,0,0,0,1,0])]
    [M3_zero, M3] = M2.children
    assert M3_zero.module == []
    assert M3_zero.prev_module == [0]
    assert M3_zero.e_images == []
    assert M3_zero.children == []
    assert M3.module == [5]
    assert M3.prev_module == [5]
    assert M3.e_images == [
        Vector([0,1,0,0,-1,0]),
    ]
    [M4] = M3.children
    assert M4 is M3

def test_exponentially_growing_Zs():
    res = ProjectiveResolution([[0,1,0,1,0,0,0],
                                [0,1,0,1,0,0,1],
                                [2,3,2,3,2,2,2],
                                [2,3,2,3,2,2,3],
                                [0,1,0,1,0,0,4],
                                [0,1,0,1,0,0,5],
                                [0,1,2,3,4,5,6]])
    assert res.homology_list(10) == [{0: 1}, {}, {0: 1}, {0: 2}, {0: 4}, {0: 8}, {0: 16}, {0: 32}, {0: 64}, {0: 128}, {0: 256}]
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    [M1] = M0.get_children()
    assert M1.module == [6]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1,-1])]
    # the kernel of *(x10-x00) on {x00,x01,x10,x11,y00,z00,1}
    # is <x00,x01,x10,x11,y00,z00>
    # Which can be covered by ZSx01 + ZSy00 + ZSz00
    [M2] = M1.children
    assert M2.module == [0, 6, 6]
    assert M2.prev_module == [6]
    assert M2.e_images == [
        Vector([0,1,0,0,0,0,0]),
        Vector([0,0,0,0,1,0,0]),
        Vector([0,0,0,0,0,1,0]),
    ]
    [M3_zero, M3] = M2.children
    assert M3_zero.module == []
    assert M3_zero.prev_module == [0]
    assert M3_zero.e_images == []
    assert M3.module == [6, 2, 2, 6, 6, 6]
    assert M3.prev_module == [6, 6]
    assert len(M3.e_images) == 6
    for row in M3.e_images:
        assert len(row) == 7 + 7
