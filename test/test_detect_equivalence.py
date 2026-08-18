from mutable_lattice import Vector
import pytest

from fast_semigroup_homology.detect_equivalence import (
    slice_resolution,
    attempt_detect_equivalence,
)

def test_slice_resolution_retract_rect22():
    op1 = [[0,1,0,1,0,0],
           [0,1,0,1,0,1],
           [2,3,2,3,2,2],
           [2,3,2,3,2,3],
           [0,1,0,1,4,4],
           [0,1,2,3,4,5]]
    op2 = [[0,1,0,1,0],
           [0,1,0,1,1],
           [2,3,2,3,2],
           [2,3,2,3,3],
           [0,1,2,3,4]]
    f = [0,1,2,3,0,4]
    res = slice_resolution(op1, op2, f)
    assert res.homology_list(10) == [{0: 1}] + [{}] * 10
    res.assert_exact()
    M0 = res.root
    assert M0.module == [5]
    assert M0.prev_module is None
    assert M0.e_images == [Vector([0,0,0,0,1])]
    assert M0.child_gen_indexes == [[0]]
    [M1] = M0.children
    assert M1.module == [4]
    assert M1.prev_module == [5]
    assert M1.e_images == [Vector([1,0,0,0,-1,0])]
    assert M1.child_gen_indexes == [[0]]
    [M2] = M1.children
    assert M2.module == [0]
    assert M2.prev_module == [4]
    assert M2.e_images == [Vector([1,0,0])]
    assert M2.child_gen_indexes == []
    assert M2.children == []

@pytest.mark.parametrize("opstring", [
    "0",
    "01;10",
    "012;120;201",
    "0123;1230;2301;3012",
    "0123;1032;2301;3210",
    "01234;12340;23401;34012;40123",
    "01010;01011;23232;23233;01234",
])
def test_detect_equivalence_id(opstring):
    op = [list(map(int, line)) for line in opstring.split(";")]
    assert attempt_detect_equivalence(op, op, list(range(len(op))))

@pytest.mark.parametrize("opstring", [
    "00;01",
    "000;012;021",
    "010;011;012",
    "000;111;012",
    "000;011;012",
])
def test_attempt_detect_equivalence_with_zero(opstring):
    op = [list(map(int, line)) for line in opstring.split(";")]
    assert attempt_detect_equivalence(op, [[0]], [0]*len(op))

def test_attempt_detect_equivalence_with_C2():
    # C2 + adjoined 1 --> C2
    assert attempt_detect_equivalence([[0,1,0],[1,0,1],[0,1,2]], [[0,1],[1,0]], [0,1,0])
    # C2 + trivially adjoined C2 --> C2
    assert attempt_detect_equivalence([[0,1,0,0],[1,0,1,1],[0,1,2,3],[0,1,3,2]], [[0,1],[1,0]], [0,1,0,0])
    # {0,1} x C2 --> C2
    assert attempt_detect_equivalence([[0,1,0,1],[1,0,1,0],[0,1,2,3],[1,0,3,2]], [[0,1],[1,0]], [0,1,0,1])
    # (left zero) x C2 + adjoined 1 --> C2
    assert attempt_detect_equivalence([[0,1,0,1,0],[1,0,1,0,1],[2,3,2,3,2],[3,2,3,2,3],[0,1,2,3,4]], [[0,1],[1,0]], [0,1,0,1,0])
    # (right zero) x C2 + adjoined 1 --> C2
    assert attempt_detect_equivalence([[0,1,2,3,0],[1,0,3,2,1],[0,1,2,3,2],[1,0,3,2,3],[0,1,2,3,4]], [[0,1],[1,0]], [0,1,0,1,0])

@pytest.mark.parametrize("opstring,f4,f5", [
    # this first map apparently isn't detected by Quillen's Thm A
    # ("0101000;0101001;2323222;2323223;0101004;0101045;0123456",0,4),
    ("0101000;0101001;2323222;2323223;0101004;0101055;0123456",4,0),
    ("0101000;0101001;2323222;2323223;0101004;2323255;0123456",4,2),
    ("0101010;0101011;2323232;2323233;0101014;2323255;0123456",4,3),
    ("0101000;0101001;2323222;2323223;0101444;0101445;0123456",0,4),
])
def test_equivalence_with_0ZZZ(opstring, f4, f5):
    op1 = [list(map(int, line)) for line in opstring.split(";")]
    op2 = [[0,1,0,1,0,0],
           [0,1,0,1,0,1],
           [2,3,2,3,2,2],
           [2,3,2,3,2,3],
           [0,1,0,1,0,4],
           [0,1,2,3,4,5]]
    f = [0,1,2,3,f4,f5,5]
    assert attempt_detect_equivalence(op1, op2, f)

@pytest.mark.parametrize("opstring,f4,f5", [
    ("0101000;0101001;2323222;2323223;0101004;0101455;0123456",0,0),
    ("0101000;0101001;2323222;2323223;0101044;0101455;0123456",0,0),
    ("0101000;0101001;2323222;2323223;0101044;2323255;0123456",0,2),
    ("0101000;0101001;2323222;2323223;0101454;0101545;0123456",0,0),
    ("0101000;0101001;2323222;2323223;0101404;0101055;0123456",0,0),
    ("0101000;0101001;2323222;2323223;0101404;2323255;0123456",0,2),
    ("0101000;0101001;2323222;2323223;0101444;0101455;0123456",0,0),
    ("0101000;0101001;2323222;2323223;0101444;0101555;0123456",0,0),
    ("0101000;0101001;2323222;2323223;0101444;2323555;0123456",0,2),
    ("0101010;0101011;2323232;2323233;0101414;2323255;0123456",0,3),
])
def test_equivalence_with_S2(opstring, f4, f5):
    op1 = [list(map(int, line)) for line in opstring.split(";")]
    op2 = [[0,1,0,1,0],
           [0,1,0,1,1],
           [2,3,2,3,2],
           [2,3,2,3,3],
           [0,1,2,3,4]]
    f = [0,1,2,3,f4,f5,4]
    assert attempt_detect_equivalence(op1, op2, f)

def test_equivalence_with_0Z202020():
    op1 = [[0,1,0,1,0,0,0],
           [0,1,0,1,0,1,1],
           [2,3,2,3,2,2,2],
           [2,3,2,3,2,3,3],
           [0,1,0,1,4,4,4],
           [0,1,2,3,4,5,6],
           [0,1,2,3,4,6,5]]
    op2 = [[0,1,0,1,0,0],
           [0,1,0,1,1,1],
           [2,3,2,3,2,2],
           [2,3,2,3,3,3],
           [0,1,2,3,4,5],
           [0,1,2,3,5,4]]
    f = [0,1,2,3,0,4,5]
    assert attempt_detect_equivalence(op1, op2, f)
