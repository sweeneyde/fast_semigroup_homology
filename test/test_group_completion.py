from fast_semigroup_homology.group_completion import (
    _minimial_ideal_structure,
    _normal_subgroup_generated,
    _quotient_group,
    group_completion,
    group_completion_functor_map,
    map_induces_isomorphism_on_group_completion,
    universal_cover_homology_list,
    pi2,
)

C1 = [[0]]

C2 = [[0,1],
      [1,0]]

C4 = [[0,1,2,3],
      [1,2,3,0],
      [2,3,0,1],
      [3,0,1,2]]

C2xC2 = [[0,1,2,3],
         [1,0,3,2],
         [2,3,0,1],
         [3,2,1,0]]

D6 = [[0,1,2,3,4,5],
      [1,2,0,5,3,4],
      [2,0,1,4,5,3],
      [3,4,5,0,1,2],
      [4,5,3,2,0,1],
      [5,3,4,1,2,0]]

rect22 = [[0,1,0,1],
          [0,1,0,1],
          [2,3,2,3],
          [2,3,2,3]]

C2_2_2_sandwich0 = [[0,1,2,3,0,1,2,3],
                    [1,0,3,2,1,0,3,2],
                    [0,1,2,3,0,1,2,3],
                    [1,0,3,2,1,0,3,2],
                    [4,5,6,7,4,5,6,7],
                    [5,4,7,6,5,4,7,6],
                    [4,5,6,7,4,5,6,7],
                    [5,4,7,6,5,4,7,6]]

C2_2_2_sandwich1 = [[0,1,2,3,0,1,2,3],
                    [1,0,3,2,1,0,3,2],
                    [0,1,2,3,1,0,3,2],
                    [1,0,3,2,0,1,2,3],
                    [4,5,6,7,4,5,6,7],
                    [5,4,7,6,5,4,7,6],
                    [4,5,6,7,5,4,7,6],
                    [5,4,7,6,4,5,6,7]]

SusC2 = [[0,1,0,1,0,1],
         [0,1,0,1,1,0],
         [2,3,2,3,2,3],
         [2,3,2,3,3,2],
         [0,1,2,3,4,5],
         [0,1,2,3,5,4]]

def adjoin_1(op):
    return [row + [i] for i, row in enumerate(op)] + [list(range(len(op)+1))]

def test_minimal_ideal_structure():
    assert _minimial_ideal_structure(C1) == ({0}, {0}, {0})
    assert _minimial_ideal_structure([[0,0],[0,1]]) == ({0}, {0}, {0})
    assert _minimial_ideal_structure([[0,1],[1,1]]) == ({1}, {1}, {1})
    assert _minimial_ideal_structure([[0,0],[1,1]]) == ({0}, {0,1}, {0})
    assert _minimial_ideal_structure([[0,1],[0,1]]) == ({0}, {0}, {0,1})
    assert _minimial_ideal_structure(rect22) == ({0}, {0,2}, {0,1})
    assert _minimial_ideal_structure(C2) == ({0,1}, {0}, {0})
    assert _minimial_ideal_structure(C2xC2) == ({0,1,2,3}, {0}, {0})
    assert _minimial_ideal_structure(C4) == ({0,1,2,3}, {0}, {0})
    assert _minimial_ideal_structure(C2_2_2_sandwich0) == ({0,1}, {0,4}, {0,2})
    assert _minimial_ideal_structure(C2_2_2_sandwich1) == ({0,1}, {0,4}, {0,2})
    assert _minimial_ideal_structure(SusC2) == ({0}, {0,2}, {0,1})

def test_normal_subgroup_generated():
    assert _normal_subgroup_generated(range(1), C1, []) == {0}
    assert _normal_subgroup_generated(range(1), C1, [0]) == {0}
    assert _normal_subgroup_generated(range(2), C2, []) == {0}
    assert _normal_subgroup_generated(range(2), C2, [0]) == {0}
    assert _normal_subgroup_generated(range(2), C2, [1]) == {0,1}
    assert _normal_subgroup_generated(range(2), C2, [0,1]) == {0,1}
    assert _normal_subgroup_generated(range(4), C4, []) == {0}
    assert _normal_subgroup_generated(range(4), C4, [0]) == {0}
    assert _normal_subgroup_generated(range(4), C4, [1]) == {0,1,2,3}
    assert _normal_subgroup_generated(range(4), C4, [2]) == {0,2}
    assert _normal_subgroup_generated(range(4), C4, [3]) == {0,1,2,3}
    assert _normal_subgroup_generated(range(4), C2xC2, []) == {0}
    assert _normal_subgroup_generated(range(4), C2xC2, [0]) == {0}
    assert _normal_subgroup_generated(range(4), C2xC2, [1]) == {0,1}
    assert _normal_subgroup_generated(range(4), C2xC2, [2]) == {0,2}
    assert _normal_subgroup_generated(range(4), C2xC2, [3]) == {0,3}
    assert _normal_subgroup_generated(range(4), C2xC2, [0,1]) == {0,1}
    assert _normal_subgroup_generated(range(4), C2xC2, [0,2]) == {0,2}
    assert _normal_subgroup_generated(range(4), C2xC2, [0,3]) == {0,3}
    assert _normal_subgroup_generated(range(4), C2xC2, [1,2]) == {0,1,2,3}
    assert _normal_subgroup_generated(range(4), C2xC2, [1,3]) == {0,1,2,3}
    assert _normal_subgroup_generated(range(4), C2xC2, [2,3]) == {0,1,2,3}
    assert _normal_subgroup_generated(range(6), D6, []) == {0}
    assert _normal_subgroup_generated(range(6), D6, [1]) == {0,1,2}
    assert _normal_subgroup_generated(range(6), D6, [2]) == {0,1,2}
    for flip in [3,4,5]:
        assert _normal_subgroup_generated(range(6), D6, [flip]) == {0,1,2,3,4,5}

def test_quotient_group():
    assert _quotient_group(range(1), C1, []) == (C1, {0:0}, [0])
    assert _quotient_group(range(2), C2, []) == (C2, {0:0, 1:1}, [0,1])
    assert _quotient_group(range(2), C2, [1]) == (C1, {0:0, 1:0}, [0])
    assert _quotient_group(range(4), C4, []) == (C4, {0:0, 1:1, 2:2, 3:3}, [0,1,2,3])
    assert _quotient_group(range(4), C4, [1]) == (C1, {0:0, 1:0, 2:0, 3:0}, [0])
    assert _quotient_group(range(4), C4, [2]) == (C2, {0:0, 1:1, 2:0, 3:1}, [0,1])
    assert _quotient_group(range(4), C2xC2, []) == (C2xC2, {0:0, 1:1, 2:2, 3:3}, [0,1,2,3])
    assert _quotient_group(range(4), C2xC2, [1]) == (C2, {0:0, 1:0, 2:1, 3:1}, [0,2])
    assert _quotient_group(range(4), C2xC2, [2]) == (C2, {0:0, 1:1, 2:0, 3:1}, [0,1])
    assert _quotient_group(range(4), C2xC2, [3]) == (C2, {0:0, 1:1, 2:1, 3:0}, [0,1])
    assert _quotient_group(range(4), C2xC2, [1,2]) == (C1, {0:0, 1:0, 2:0, 3:0}, [0])
    assert _quotient_group(range(6), D6, []) == (D6, {0:0, 1:1, 2:2, 3:3, 4:4, 5:5}, [0,1,2,3,4,5])
    assert _quotient_group(range(6), D6, [1]) == (C2, {0:0, 1:0, 2:0, 3:1, 4:1, 5:1}, [0,3])
    assert _quotient_group(range(6), D6, [3]) == (C1, {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}, [0])

def test_group_completion():
    assert group_completion(C1) == (C1, [0], [0])
    assert group_completion(C2) == (C2, [0,1], [0,1])
    assert group_completion(C2xC2) == (C2xC2, [0,1,2,3], [0,1,2,3])
    assert group_completion(C4) == (C4, [0,1,2,3], [0,1,2,3])
    assert group_completion(D6) == (D6, [0,1,2,3,4,5], [0,1,2,3,4,5])
    assert group_completion([[0,0],[1,1]]) == (C1, [0,0], [0])
    assert group_completion([[0,1],[0,1]]) == (C1, [0,0], [0])
    assert group_completion(C2_2_2_sandwich0) == (C2, [0,1,0,1,0,1,0,1], [0,1])
    assert group_completion(C2_2_2_sandwich1) == (C1, [0,0,0,0,0,0,0,0], [0])

def test_group_completion_functor_map():
    assert group_completion_functor_map(C1, C1, [0]) == (C1, C1, [0])
    assert group_completion_functor_map(C2, C1, [0,0]) == (C2, C1, [0,0])
    assert group_completion_functor_map(C1, C2, [0]) == (C1, C2, [0])
    assert group_completion_functor_map(C2, C2, [0,0]) == (C2, C2, [0,0])
    assert group_completion_functor_map(C2, C2, [0,1]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2, C4, [0,2]) == (C2, C4, [0,2])
    assert group_completion_functor_map(C2, C4, [0,0]) == (C2, C4, [0,0])
    assert group_completion_functor_map(C4, C2, [0,1,0,1]) == (C4, C2, [0,1,0,1])
    assert group_completion_functor_map(C4, C2, [0,1,0,1]) == (C4, C2, [0,1,0,1])
    assert group_completion_functor_map(C2, C2_2_2_sandwich0, [0,1]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2, C2_2_2_sandwich0, [2,3]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2, [0,1,0,1,0,1,0,1]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2, [0,0,0,0,0,0,0,0]) == (C2, C2, [0,0])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,2,3,4,5,6,7]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,0,2,2,4,4,6,6]) == (C2, C2, [0,0])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,2,3,0,1,2,3]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,0,1,4,5,4,5]) == (C2, C2, [0,1])
    assert group_completion_functor_map(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,0,1,0,1,0,1]) == (C2, C2, [0,1])

def test_map_induces_isomorphism_on_group_completion():
    assert map_induces_isomorphism_on_group_completion(C1, C1, [0])
    assert not map_induces_isomorphism_on_group_completion(C2, C1, [0,0])
    assert not map_induces_isomorphism_on_group_completion(C1, C2, [0])
    assert not map_induces_isomorphism_on_group_completion(C2, C2, [0,0])
    assert map_induces_isomorphism_on_group_completion(C2, C2, [0,1])
    assert map_induces_isomorphism_on_group_completion(C2, C2_2_2_sandwich0, [0,1])
    assert map_induces_isomorphism_on_group_completion(C2, C2_2_2_sandwich0, [2,3])
    assert map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2, [0,1,0,1,0,1,0,1])
    assert not map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2, [0,0,0,0,0,0,0,0])
    assert map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,2,3,4,5,6,7])
    assert not map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,0,2,2,4,4,6,6])
    assert map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,2,3,0,1,2,3])
    assert map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,0,1,4,5,4,5])
    assert map_induces_isomorphism_on_group_completion(C2_2_2_sandwich0, C2_2_2_sandwich0, [0,1,0,1,0,1,0,1])

def test_universal_cover_homology_list():
    assert universal_cover_homology_list(C1, 6) == [{0: 1}, {}, {}, {}, {}, {}, {}]
    assert universal_cover_homology_list(C2, 6) == [{0: 1}, {}, {}, {}, {}, {}, {}]
    assert universal_cover_homology_list(C2xC2, 6) == [{0: 1}, {}, {}, {}, {}, {}, {}]
    assert universal_cover_homology_list(D6, 6) == [{0: 1}, {}, {}, {}, {}, {}, {}]
    assert universal_cover_homology_list(adjoin_1(rect22), 6) == [{0: 1}, {}, {0: 1}, {}, {}, {}, {}]
    # The universal cover of B(C2)vS2 is E(C2)vS2vS2 (one S2 at each vertex of E(C2))
    assert universal_cover_homology_list(adjoin_1(C2_2_2_sandwich0), 6) == [{0: 1}, {}, {0: 2}, {}, {}, {}, {}]

def test_pi2():
    assert pi2(C1) == {}
    assert pi2(C2) == {}
    assert pi2(C2xC2) == {}
    assert pi2(D6) == {}
    assert pi2(adjoin_1(rect22)) == {0: 1}
    assert pi2(adjoin_1(C2_2_2_sandwich0)) == {0: 2}
    assert pi2(adjoin_1(C2_2_2_sandwich1)) == {0: 1}
    assert pi2(SusC2) == {2: 1}
