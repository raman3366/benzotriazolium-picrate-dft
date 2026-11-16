"""Write a CSV copy of every table in tables/ (tbl*.md -> tbl*.csv).

    python scripts/tables/tables_to_csv.py

The Markdown files (pandoc pipe tables, made by make_tables.py) are the record; the CSV copies are for reading
the numbers by program. Each CSV has the header and every row of the table in order, with Markdown emphasis and
code marks removed, the minus sign U+2212 written as '-', and no-break and thin spaces written as ordinary
spaces. Rows that only name a block of the table (for example "Bond lengths (Å), picrate") keep their text in the
first column and leave the others empty. Sub- and superscript marks (~ ~, ^ ^) are dropped. The caption, with the notes, stays in the Markdown file.
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths  # noqa: E402


def cells(line):
    parts = re.split(r'(?<!\\)\|', line.strip()[1:-1])           # pipes escaped as \| belong to the cell
    return [clean(p.replace('\\|', '|')) for p in parts]


def clean(s):
    s = s.strip().replace('−', '-').replace(' ', ' ').replace(' ', ' ').replace(' ', ' ')
    s = re.sub(r'\*\*([^*]+?)\*\*', r'\1', s)                   # **strong**
    s = re.sub(r'\*([^*]+?)\*', r'\1', s)                         # *emphasis*
    s = re.sub(r'~([^~]+?)~', r'\1', s)                           # ~subscript~
    s = re.sub(r'\^([^^]+?)\^', r'\1', s)                        # ^superscript^
    s = s.replace('`', '')
    return s


def convert(md):
    lines = md.read_text(encoding='utf-8').splitlines()
    rows = [l for l in lines if l.startswith('|')]
    if len(rows) < 2 or not re.fullmatch(r'\|[\s:|-]+\|', rows[1].strip()):
        raise SystemExit(f'{md.name}: no pipe table found')
    out = [cells(rows[0])] + [cells(r) for r in rows[2:]]
    n = len(out[0])
    for r in out:
        if len(r) != n:
            raise SystemExit(f'{md.name}: row with {len(r)} cells, header has {n}: {r}')
    dst = md.with_suffix('.csv')
    with open(dst, 'w', newline='', encoding='utf-8') as f:
        csv.writer(f, lineterminator='\n').writerows(out)
    return dst, len(out) - 1, n


if __name__ == '__main__':
    for md in sorted(paths.TABLES.glob('tbl*.md')):
        dst, nr, nc = convert(md)
        print(f'{dst.relative_to(paths.REPO)}: {nr} rows x {nc} columns')
