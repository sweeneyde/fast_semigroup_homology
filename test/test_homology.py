import fast_monoid_homology.kernels
import doctest

def test_homology_doctest():
    results = doctest.testmod(fast_monoid_homology.kernels)
    assert results.failed == 0
