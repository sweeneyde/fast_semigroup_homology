import fast_semigroup_homology.kernels
import doctest

def test_homology_doctest():
    results = doctest.testmod(fast_semigroup_homology.kernels)
    assert results.failed == 0
