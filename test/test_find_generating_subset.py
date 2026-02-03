from fast_monoid_homology.find_generating_subset import (
    find_generating_subset
)
from mutable_lattice import Vector

def test_find_generating_subset():
    # cover the kernel of
    # Z <---augmentation--- ZM
    # where M is 2x2 rectangular band, plus unit
    assert find_generating_subset(
        5,
        [
            Vector([1,0,0,0,-1]), # x00 - 1
            Vector([0,1,0,0,-1]), # x01 - 1
            Vector([0,0,1,0,-1]), # x10 - 1
            Vector([0,0,0,1,-1]), # x11 - 1
        ],
        [
            Vector([0,1,0,1,0]),
            Vector([0,1,0,1,1]),
            Vector([2,3,2,3,2]),
            Vector([2,3,2,3,3]),
            Vector([0,1,2,3,4]),
        ],
        [5, 5, 5, 5],
        True,
        True,
    ) == [
        Vector([1,0,0,0,-1]), # x00 - 1
        Vector([0,0,1,0,-1]), # x10 - 1
    ]

    # cover the kernel of
    # Z <---augmentation--- ZM*x00
    assert find_generating_subset(
        2,
        [
            Vector([1,-1]), # x00-x10
        ],
        [
            Vector([0,0]),
            Vector([0,0]),
            Vector([1,1]),
            Vector([1,1]),
            Vector([0,1]),
        ],
        [5],
        True,
        False,
    ) == [
        Vector([1,-1]), # x00-x10
    ]

    assert find_generating_subset(
        1,
        [Vector([1])],
        [Vector([0]), Vector([0]), Vector([0]), Vector([0]), Vector([0])],
        [1], False, False,
    ) == [
        Vector([1])
    ]
