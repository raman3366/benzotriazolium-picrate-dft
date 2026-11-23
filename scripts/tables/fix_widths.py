"""Set the separator-row dash counts of the revision tables from their content, so that pandoc's relative column
widths give numbers room and let long headers wrap (pandoc takes column widths from the dashes when a row is wider
than --columns). Run after the rev_*_tables.py scripts:  python scripts/tables/fix_widths.py"""
import re, sys, unicodedata
from pathlib import Path

def vis(s):
    s = re.sub(r'[*~^]', '', s)
    return sum(0 if unicodedata.combining(c) else 1 for c in s.strip())

def fix(p):
    lines = p.read_text().split('\n')
    rows = [i for i, l in enumerate(lines) if l.startswith('|')]
    if len(rows) < 3:
        return False
    sep = rows[1]
    cells = lambda l: [c for c in l.strip().strip('|').split('|')]
    head = cells(lines[rows[0]])
    body = [cells(lines[i]) for i in rows[2:]]
    seps = cells(lines[sep])
    widths = []
    numre = re.compile(r'^[\s−+\-–.,0-9()%^a-z]*$')
    for j in range(len(seps)):
        # group rows (text in the first cell only) may wrap; they do not set the width
        col = [r[j] for r in body if j < len(r) and any(c.strip() for c in r[1:])]
        cmax = max((vis(c) for c in col), default=1)
        hword = max((vis(w) for w in re.split(r'\s+', head[j]) if w), default=1)
        textcol = sum(1 for c in col if c.strip() and not numre.match(c.strip())) > len(col) / 2
        w = max(min(cmax, 22) if textcol else cmax, hword, 3) + 1
        widths.append(w)
    new = []
    for j, s in enumerate(seps):
        s = s.strip(); left, right = s.startswith(':'), s.endswith(':')
        n = widths[j] - left - right
        new.append((':' if left else '') + '-' * n + (':' if right else ''))
    lines[sep] = '|' + '|'.join(new) + '|'
    p.write_text('\n'.join(lines))
    return True

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths                                                     # noqa: E402
here = paths.TABLES                                              # dataset: tables/
for p in sorted(here.glob('tbl*_rev_*.md')):
    print(p.name, fix(p))
