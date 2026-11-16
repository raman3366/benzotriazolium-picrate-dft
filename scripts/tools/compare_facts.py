"""Compare two facts.json files key by key.

    python scripts/tools/compare_facts.py <reference facts.json> <new facts.json>

Prints the number of keys, the keys whose value differs (value and type are compared exactly) and, for
information, the keys whose unit, case, source or note differs. Exit status 1 if any value differs or a key is
missing on either side.
"""
import json
import sys


def main(a_path, b_path):
    a = json.load(open(a_path, encoding='utf-8'))
    b = json.load(open(b_path, encoding='utf-8'))
    only_a = [k for k in a if k not in b]
    only_b = [k for k in b if k not in a]
    common = [k for k in a if k in b]
    dv = [k for k in common if type(a[k]['value']) is not type(b[k]['value']) or a[k]['value'] != b[k]['value']]
    other = {}
    for k in common:
        for f in sorted(set(a[k]) | set(b[k])):
            if f != 'value' and a[k].get(f) != b[k].get(f):
                other.setdefault(f, []).append(k)
    print(f'{len(a)} / {len(b)} keys; {len(dv)} values differ; {len(only_a)} only in the first, {len(only_b)} only in the second')
    for k in dv:
        print(f'  value differs: {k}: {a[k]["value"]!r} -> {b[k]["value"]!r}')
    for k in only_a + only_b:
        print(f'  missing on one side: {k}')
    for f, ks in other.items():
        print(f'  {f} differs in {len(ks)} keys: {", ".join(ks[:8])}{" ..." if len(ks) > 8 else ""}')
    return 1 if (dv or only_a or only_b) else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2]))
