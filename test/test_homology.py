import fast_semigroup_homology.homology
import doctest

def test_homology_doctest():
    results = doctest.testmod(fast_semigroup_homology.homology)
    assert results.attempted > 0
    assert results.failed == 0
