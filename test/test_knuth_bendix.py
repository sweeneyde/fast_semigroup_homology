from fast_semigroup_homology.cup_products.rewriting.knuth_bendix import (
    shortlex_ordered,
    normalize_rules,
    reduced,
    get_critical_pairs_from_normalized,
    kb_complete,
    kb_normalize,
)

def test_doctest_complete():
    rules = [("xxx", ""), ("yyy", ""), ("xyxyxy", "")]
    assert normalize_rules(rules) == rules
    criticals = get_critical_pairs_from_normalized(rules)
    assert criticals == [("xyxyx", "yy"), ("yxyxy", "xx")]
    rules2 = normalize_rules(rules + criticals)
    assert rules2 == [
        ("xxx", ""), ("yyy", ""), ("xyxyx", "yy"), ("yxyxy", "xx")
    ]
    criticals2 = get_critical_pairs_from_normalized(rules2)
    assert criticals2 == [("yxyx", "xxyy"), ("yyxx", "xyxy")]
    rules3 = normalize_rules(rules2 + criticals2)
    assert rules3 == [
        ("xxx", ""), ("yyy", ""), ("yxyx", "xxyy"), ("yyxx", "xyxy")
    ]
    assert normalize_rules(rules3) == rules3
    assert get_critical_pairs_from_normalized(rules3) == []
    assert kb_complete(rules) == rules3
    assert rules3 == [('xxx', ''), ('yyy', ''), ('yxyx', 'xxyy'), ('yyxx', 'xyxy')]
