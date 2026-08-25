"""
Apply the Knuth-Bendix algorithm for monoid presentations.

Given an input list of pair of strings that are declared equal,
produce an output list of such pairs of strings so that
two strings are equal according to the input list
iff they are equal according to the output list.

Example:

    >>> kb_normalize('xy', [('xxx', ''), ('yyy', ''), ('xyxyxy', '')])
    [('xxx', ''), ('yyy', ''), ('yxyx', 'xxyy'), ('yyxx', 'xyxy')]

The output list is normalized so that applying
`word = word.replace(left, right)` for each `(left, right)` pair
repeatedly in any order will eventually terminate, and will
always reach the same fully-reduced final answer regardless of the
order in which the rules are applied.
Thus, the word problem for the monoid is solved: two strings
are equal if and only if they reduce to the same word.

The input list in the example is not normalized because
starting with xxxyxyxy, we can use the input list of replacements
to produce yxyxy, which is then irreducible by the input list,
while we could also have used a different rule to produce xx,
which is likewise irreducible; order mattered.
"""

__all__ = ["normalize"]

from collections import defaultdict, deque

def shortlex_ordered(a: str, b: str):
    """Return a pair of strings sorted into shortlex order"""
    if (len(a), a) > (len(b), b):
        return (a, b)
    else:
        return (b, a)

def normalize_rules(rules):
    # Ensure that no right-hand sides contain left-hand sides.
    # Also ensure that no left-hand sides contain other left-hand sides.
    normal = []
    unnormal = deque(rules)
    while unnormal:
        next_left, next_right = shortlex_ordered(*unnormal.popleft())
        if next_left == next_right:
            continue
        assert next_left not in next_right # because of shortlex
        if (next_left, next_right) in normal:
            continue
        for left, right in normal:
            if left in next_left:
                unnormal.appendleft((next_left.replace(left, right), next_right))
                break
            if left in next_right:
                unnormal.appendleft((next_left, next_right.replace(left, right)))
                break
            if next_left in left:
                normal.remove((left, right))
                unnormal.append((left.replace(next_left, next_right), right))
                unnormal.appendleft((next_left, next_right))
                break
            if next_left in right:
                normal.remove((left, right))
                unnormal.append((left, right.replace(next_left, next_right)))
                unnormal.appendleft((next_left, next_right))
                break
        else:
            # no break: this rule passed all checks for being normal.
            normal.append((next_left, next_right))
    return normal

def reduced(word, rule_list):
    # This is only called when rule_list has shortlex-shrinking words
    while True:
        word0 = word
        for left, right in rule_list:
            word = word.replace(left, right)
        if word == word0:
            return word

def get_critical_pairs_from_normalized(rules):
    # Phase 2: finding critical pairs
    # Since we've normalized the rules, only consider
    # (AB-->..., BC-->...) cases, not (ABC-->...,B-->...) cases.
    prefix_to_rules = defaultdict(list)
    suffix_to_rules = defaultdict(list)
    for left, right in rules:
        for i in range(1, len(left)):
            prefix, suffix = left[:i], left[i:]
            prefix_to_rules[prefix].append((suffix, right))
            suffix_to_rules[suffix].append((prefix, right))
    pairs = set()
    for middle in prefix_to_rules.keys() & suffix_to_rules.keys():
        prefix_rules = prefix_to_rules[middle]
        suffix_rules = suffix_to_rules[middle]
        for suffix_rule_prefix, suffix_rule_right in suffix_rules:
            for prefix_rule_suffix, prefix_rule_right in prefix_rules:
                result1 = reduced(suffix_rule_right + prefix_rule_suffix, rules)
                result2 = reduced(suffix_rule_prefix + prefix_rule_right, rules)
                if result1 != result2:
                    pairs.add(shortlex_ordered(result1, result2))
    return sorted(pairs)

def kb_complete(rules, *, iteration_limit=20):
    # copy to verify at the end
    rules0 = [(left, right) for left, right in rules]
    iterations = 0
    while True:
        rules = normalize_rules(rules)
        critical_pairs = list(get_critical_pairs_from_normalized(rules))
        if not critical_pairs:
            break
        rules.extend(critical_pairs)
        iterations += 1
        if iterations >= iteration_limit:
            raise RuntimeError(f"Did not normalize after {iteration_limit} loops")
    for left0, right0 in rules0:
        assert reduced(left0, rules) == reduced(right0, rules)
    return rules

def kb_normalize(alphabet, rules0, **kwargs):
    rules = kb_complete(rules0, **kwargs)
    # remove single-letter left-hand sides.
    normal_rules = []
    for rule in rules:
        if len(rule[0]) == 1:
            alphabet = alphabet.replace(rule[0][0], '')
        else:
            normal_rules.append(rule)
    return (alphabet, normal_rules)
