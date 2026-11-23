"""Revision facts: the experimental section (co-authors' values, R1.2) and quantities derived from them."""
import json, re
from .common import *

# standard atomic weights (IUPAC 2021, abridged), used only to convert the co-authors' masses to amounts
AW = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999}
def mr(formula):
    return sum(AW[e] * (int(n) if n else 1) for e, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula))

def build(add, facts):
    src = 'papers/btpc/data/experimental.json'
    d = json.load(open(ROOT / src))
    why = 'co-authors\' answer, item {} (answers to the authors\' questionnaire of the revision, not deposited)'
    for k, v in d.items():
        if k.startswith('_'):
            continue
        add(f'rev.exp.{k}', v['value'], v['unit'], None, f'{src}:{k}', why.format(v['item']) + ('; ' + v['note'] if v.get('note') else ''))
    m_bta, m_pa = mr('C6H5N3'), mr('C6H3N3O7')
    n_bta = d['benzotriazole_mass_g']['value'] / m_bta * 1000
    n_pa = d['picric_acid_mass_g']['value'] / m_pa * 1000
    add('rev.exp.benzotriazole_mmol', n_bta, 'mmol', None, src, f'derived: 5.00 g / M(C6H5N3) = {m_bta:.3f} g mol-1 (IUPAC atomic weights)')
    add('rev.exp.picric_acid_mmol_max', n_pa, 'mmol', None, src,
        f'derived: 3.86 g / M(C6H3N3O7) = {m_pa:.3f} g mol-1; an upper bound, because the acid was weighed moist')
    add('rev.exp.molar_ratio_min', n_bta / n_pa, '', None, src, 'derived: benzotriazole : picric acid, lower bound (moist acid)')
