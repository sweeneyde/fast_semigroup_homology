from fast_semigroup_homology.projective_resolution import ProjectiveResolution
from mutable_lattice import Vector

def test_trivial():
    res = ProjectiveResolution([[0]])
    assert res.homology_list(10) == [{0: 1}] + [{}] * 10
    res.assert_exact()
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    assert M0.child_gen_indexes == []
    assert M0.children == []

def test_C2():
    res = ProjectiveResolution([[0,1],[1,0]])
    assert res.homology_list(10) == [{0: 1}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}]
    res.assert_exact()
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    assert M0.child_gen_indexes == [[0]]
    [M1] = M0.children
    assert M1.module == [0]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1, -1])]
    assert M1.child_gen_indexes == [[0]]
    [M2] = M1.children
    assert M2.module == [0]
    assert M2.prev_module == [0]
    assert M2.e_images == [Vector([1, 1])]
    assert M2.child_gen_indexes == [[0]]
    [M3] = M2.children
    assert M3 is M1

def test_C3():
    res = ProjectiveResolution([[0,1,2],[1,2,0],[2,0,1]])
    assert res.homology_list(10) == [{0: 1}] + [{3: 1}, {}] * 5
    res.assert_exact()
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    assert M0.child_gen_indexes == [[0]]
    [M1] = M0.children
    assert M1.module == [0]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1,0,-1])]
    assert M1.child_gen_indexes == [[0]]
    [M2] = M1.children
    assert M2.module == [0]
    assert M2.prev_module == [0]
    assert M2.e_images == [Vector([1, 1, 1])]
    assert M2.child_gen_indexes == [[0]]
    [M3] = M2.children
    assert M3 is M1

def test_rect22():
    res = ProjectiveResolution([[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]])
    assert res.homology_list(10) == [{0: 1}, {}, {0: 1}, {}, {}, {}, {}, {}, {}, {}, {}]
    res.assert_exact()
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    assert M0.child_gen_indexes == [[0]]
    [M1] = M0.children
    assert M1.module == [4]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1, -1])]
    assert M1.child_gen_indexes == [[0]]
    [M2] = M1.children
    assert M2.module == [0, 0]
    assert M2.prev_module == [4]
    assert M2.e_images == [Vector([1,0,0,0,0]),
                           Vector([0,1,0,0,0])]
    assert M2.child_gen_indexes == []
    assert M2.children == []

def test_infinitely_many_Zs():
    res = ProjectiveResolution([[0,1,0,1,0,0],[0,1,0,1,0,1],[2,3,2,3,2,2],[2,3,2,3,2,3],[0,1,0,1,0,4],[0,1,2,3,4,5]])
    assert res.homology_list(10) == [{0: 1}, {}] + [{0: 1}] * 9
    res.assert_exact()
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    assert M0.child_gen_indexes == [[0]]
    [M1] = M0.children
    assert M1.module == [5]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1,-1])]
    assert M1.child_gen_indexes == [[0]]
    [M2] = M1.children
    assert M2.module == [0, 5]
    assert M2.prev_module == [5]
    assert M2.e_images == [Vector([0,1,0,0,0,0]),
                           Vector([0,0,0,0,1,0])]
    assert M2.child_gen_indexes == [[1]]
    [M3] = M2.children
    assert M3.module == [5]
    assert M3.prev_module == [5]
    assert M3.e_images == [
        Vector([0,1,0,0,-1,0]),
    ]
    assert M3.child_gen_indexes == [[0]]
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
    # The following would not be fast enough if we didn't use homology_with_shift:
    assert res.homology_list(102)[-1] == {0: 2**100}
    assert res.homology_list(1002)[-1] == {0: 2**1000}
    res.assert_exact()
    M0 = res.root
    assert M0.module == [0]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([1])]
    assert M0.child_gen_indexes == [[0]]
    [M1] = M0.children
    assert M1.module == [6]
    assert M1.prev_module == [0]
    assert M1.e_images == [Vector([1,-1])]
    assert M1.child_gen_indexes == [[0]]
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
    assert M2.child_gen_indexes == [[1, 2]]
    # The map ZSx01(+)ZS(+)ZS ---(*x01,*y00,*z00)---> ZS sends:
    #     (x01,0,0) to x01
    #     (x11,0,0) to x11
    #     (0,1,0) to y00
    #     (0,x00,0),(0,x01,0),(0,y00,0),(0,z00,0) to x00
    #     (0,x10,0),(0,x11,0) to x10
    #     (0,0,1) to z00
    #     (0,0,x00),(0,0,x01),(0,0,y00),(0,0,z00) to x00
    #     (0,0,x10),(0,0,x11) to x10
    # The kernel of this can be covered in various ways. Here's one solution:
    [M3] = M2.children
    assert M3.module == [6, 2, 6, 6, 6]
    assert M3.prev_module == [6, 6] # The kernel is trivial on the first summand.
    assert M3.e_images == [
        Vector([0,1,0,0,0,0,0, 0,0,0,0,0,-1,0]),
        Vector([0,0,1,0,0,0,0, 0,0,0,-1,0,0,0]),
        Vector([0,0,0,0,1,0,0, 0,0,0,0,0,-1,0]),
        Vector([0,0,0,0,0,1,0, 0,0,0,0,0,-1,0]),
        Vector([0,0,0,0,0,0,0, 0,0,0,0,1,-1,0]),
    ]
    assert M3.child_gen_indexes == [[0], [2, 3], [4]]
    [M4a, M4b, M4c] = M3.children
    assert M4b is M3
    assert M4c is M2
    assert M4a.module == [6, 6]
    assert M4a.prev_module == [6]
    assert M4a.e_images == [
        Vector([0,1,0,0,0,-1,0]),
        Vector([0,0,0,0,1,-1,0]),
    ]
    assert M4a.child_gen_indexes == [[0], [1]]
    [M5a, M5b] = M4a.children
    assert M5a is M4a
    assert M5b is M2

def test_suspended_C2():
    res = ProjectiveResolution([[0,1,0,1,0,0],
                                [0,1,0,1,1,1],
                                [2,3,2,3,2,2],
                                [2,3,2,3,3,3],
                                [0,1,2,3,4,5],
                                [2,3,0,1,5,4]])
    assert res.homology_list(10) == [{0: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}, {}, {2: 1}]
    assert res.homology_list(100) == [{0: 1}] + [{}, {2: 1}] * 50
    assert res.homology_list(1000) == [{0: 1}] + [{}, {2: 1}] * 500
    res.assert_exact()

def test_C2xC2():
    res = ProjectiveResolution([[0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]])
    assert res.homology_list(10) == [{0: 1}, {2: 2}, {2: 1}, {2: 3}, {2: 2}, {2: 4}, {2: 3}, {2: 5}, {2: 4}, {2: 6}, {2: 5}]
    res.assert_exact()

def test_S3():
    res = ProjectiveResolution([[0,1,2,3,4,5],[1,0,3,2,5,4],[2,5,4,1,0,3],[3,4,5,0,1,2],[4,3,0,5,2,1],[5,2,1,4,3,0]])
    assert res.homology_list(10) == [{0: 1}, {2: 1}, {}, {6: 1}, {}, {2: 1}, {}, {6: 1}, {}, {2: 1}, {}]
    res.assert_exact()

def test_D8():
    op = [[0,1,2,3,4,5,6,7],[1,0,4,5,2,3,7,6],[2,7,0,6,5,4,3,1],[3,5,6,0,7,1,2,4],[4,6,1,7,3,2,5,0],[5,3,7,1,6,0,4,2],[6,4,3,2,1,7,0,5],[7,2,5,4,0,6,1,3]]
    res = ProjectiveResolution(op)
    assert res.homology_list(10) == [{0: 1}, {2: 2}, {2: 1}, {2: 2, 4: 1}, {2: 2}, {2: 4}, {2: 3}, {2: 4, 4: 1}, {2: 4}, {2: 6}, {2: 5}]
    res.assert_exact()

def test_Q8():
    op = [[0,1,2,3,4,5,6,7],[1,3,4,5,6,0,7,2],[2,7,3,6,1,4,0,5],[3,5,6,0,7,1,2,4],[4,2,5,7,3,6,1,0],[5,0,7,1,2,3,4,6,],[6,4,0,2,5,7,3,1],[7,6,1,4,0,2,5,3]]
    res = ProjectiveResolution(op)
    assert res.homology_list(10) == [{0: 1}, {2: 2}, {}, {8: 1}, {}, {2: 2}, {}, {8: 1}, {}, {2: 2}, {}]
    res.assert_exact()

def test_huge_torsion():
    op = [[0,1,0,1,0,0,0,0,0,0,0],
          [0,1,0,1,0,0,0,0,0,0,1],
          [2,3,2,3,2,2,2,2,2,2,2],
          [2,3,2,3,2,2,2,2,2,2,3],
          [0,1,0,1,0,0,0,0,0,0,4],
          [0,1,0,1,0,0,0,0,0,0,5],
          [0,1,0,1,0,0,0,0,0,0,6],
          [0,1,0,1,0,0,0,4,4,5,7],
          [0,1,0,1,0,0,0,6,5,6,8],
          [0,1,0,1,0,0,0,0,6,4,9],
          [0,1,2,3,4,5,6,7,8,9,10]]
    res = ProjectiveResolution(op)
    assert res.homology_list(6) == [{0: 1},
                                    {},
                                    {0: 1},
                                    {0: 3},
                                    {0: 6},
                                    {0: 9},
                                    {0: 9, 1494640: 1}]
    res.assert_exact()

def test_C2_modules():
    op = [[0, 1], [1, 0]]
    tor_Z_Z = ProjectiveResolution(op).homology_list(10)
    assert tor_Z_Z == [{0: 1}] + [{2: 1}, {}] * 5
    tor_ZC2_Z = ProjectiveResolution(op).homology_list(10, right_S_set_action=op)
    assert tor_ZC2_Z == [{0: 1}] + [{}] * 10
    tor_Z_ZC2 = ProjectiveResolution(op, left_S_set_action=op).homology_list(10)
    assert tor_Z_ZC2 == [{0: 1}] + [{}] * 10
    tor_Z_ZZ = ProjectiveResolution(op, left_S_set_action=[[0,1],[0,1]]).homology_list(10)
    assert tor_Z_ZZ == [{0: 2}] + [{2: 2}, {}] * 5
    tor_ZZ_Z = ProjectiveResolution(op).homology_list(10, right_S_set_action=[[0,0],[1,1]])
    assert tor_ZZ_Z == [{0: 2}] + [{2: 2}, {}] * 5
    tor_ZZ_ZZ = ProjectiveResolution(op, left_S_set_action=[[0,1],[0,1]]).homology_list(10, right_S_set_action=[[0,0],[1,1]])
    assert tor_ZZ_ZZ == [{0: 4}] + [{2: 4}, {}] * 5
    tor_0_Z = ProjectiveResolution(op).homology_list(10, right_S_set_action=[])
    assert tor_0_Z == [{}] * 11
    tor_Z_0 = ProjectiveResolution(op, left_S_set_action=[[], []]).homology_list(10)
    assert tor_Z_0 == [{}] * 11

def test_C4_modules():
    op = [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]
    tor_Z_Z = ProjectiveResolution(op).homology_list(10)
    assert tor_Z_Z == [{0: 1}] + [{4: 1}, {}] * 5
    tor_ZC4_Z = ProjectiveResolution(op).homology_list(10, right_S_set_action=op)
    assert tor_ZC4_Z == [{0: 1}] + [{}] * 10
    tor_ZC2_Z = ProjectiveResolution(op).homology_list(10, right_S_set_action=[[0, 1, 0, 1], [1, 0, 1, 0]])
    assert tor_ZC2_Z == [{0: 1}] + [{2: 1}, {}] * 5
    tor_Z_ZC4 = ProjectiveResolution(op, left_S_set_action=op).homology_list(10)
    assert tor_Z_ZC4 == [{0: 1}] + [{}] * 10
    tor_Z_ZC2 = ProjectiveResolution(op, left_S_set_action=[[0,1],[1,0],[0,1],[1,0]]).homology_list(10)
    assert tor_Z_ZC2 == [{0: 1}] + [{2: 1}, {}] * 5
    tor_ZC2_ZC2 = ProjectiveResolution(op, left_S_set_action=[[0,1],[1,0],[0,1],[1,0]]).homology_list(10, right_S_set_action=[[0,1,0,1],[1,0,1,0]])
    assert tor_ZC2_ZC2 == [{0: 2}] + [{2: 2}, {}] * 5

def test_rect22_modules():
    op = [[0,1,0,1,0],[0,1,0,1,1],[2,3,2,3,2],[2,3,2,3,3],[0,1,2,3,4]]
    tor_Z_Z = ProjectiveResolution(op).homology_list(10)
    assert tor_Z_Z == [{0: 1}, {}, {0: 1}, {}, {}, {}, {}, {}, {}, {}, {}]
    tor_eZS_Z = ProjectiveResolution(op).homology_list(10, right_S_set_action=[[0,1,0,1,0],[0,1,0,1,1]])
    assert tor_eZS_Z == [{0: 1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
    tor_Z_ZSe = ProjectiveResolution(op, left_S_set_action=[[0,0],[0,0],[1,1],[1,1],[0,1]]).homology_list(10)
    assert tor_Z_ZSe == [{0: 1}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}]
