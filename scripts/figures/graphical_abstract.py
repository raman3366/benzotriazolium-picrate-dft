"""Graphical abstract (a data plot drawn by this script from the computed data).

Left: gas-phase proton potential along the relaxed B3LYP scan, B3LYP and DF-MP2, with the B3LYP proton ground-state
density (data/derived/rev_proton_dvr_grid.csv). Right: IEF-PCM (water) E(neutral) - E(salt) by method
(data/derived/rev_proton_pcm.csv and facts.json rev.proton.cc.best_kcal / best_range_kcal).
Output: figures/graphical_abstract/graphical_abstract.{png,pdf,tif}, 13 x 5 cm at 600 dpi (>= 1328 x 531 px).

    python scripts/figures/graphical_abstract.py      (from the repository root; writes into figures/graphical_abstract/)
"""
import csv, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = paths.legacy_root(quiet=True) / 'papers' / 'btpc'          # dataset: the view of scripts/paths.py
D = P / 'data' / 'derived'
F = json.load(open(P / 'data' / 'facts.json'))
OUT = P / 'figures' / 'graphical_abstract'
OUT.mkdir(exist_ok=True)

grid = list(csv.DictReader(open(D / 'rev_proton_dvr_grid.csv')))
pcm = {r['method']: float(r['dE_neutral_minus_salt_kcal']) for r in csv.DictReader(open(D / 'rev_proton_pcm.csv'))}
cc = F['rev.proton.cc.best_kcal']['value']; lo, hi = F['rev.proton.cc.best_range_kcal']['value']
assert abs(pcm['b3lyp'] - F['rev.proton.pcm.b3lyp.dE_neutral_minus_salt_kcal']['value']) < 1e-3

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7})
fig = plt.figure(figsize=(13 / 2.54, 5 / 2.54))
ax = fig.add_axes([0.075, 0.2, 0.36, 0.62]); bx = fig.add_axes([0.64, 0.2, 0.335, 0.62])
r = [float(g['r_HO_A']) for g in grid if g['in_scan'] == '1']
for key, col, lab in (('b3lyp', '#d97706', 'B3LYP'), ('dfmp2', '#7c3aed', 'DF-MP2')):
    v = [float(g['V_kcal_' + key]) for g in grid if g['in_scan'] == '1']
    ax.plot(r, v, color=col, lw=1.4, label=lab)
psi = [float(g['psi2_v0_b3lyp']) for g in grid if g['in_scan'] == '1']
m = max(psi); ax.fill_between(r, 0, [2.2 * p / m for p in psi], color='#d97706', alpha=0.25, lw=0)
ax.set_xlim(0.9, 1.7); ax.set_ylim(-0.3, 6); ax.invert_xaxis()
ax.set_xlabel('r(H···O1) (Å)'); ax.set_ylabel('E (kcal mol$^{-1}$)')
ax.legend(frameon=False, loc='center left', bbox_to_anchor=(0.0, 0.62), fontsize=6.5)
ax.text(0.03, 0.95, 'H on N4 (salt)', transform=ax.transAxes, ha='left', va='top', fontsize=6)
ax.text(0.97, 0.95, 'H on O1\n(neutral)', transform=ax.transAxes, ha='right', va='top', fontsize=6)
ax.set_title('Gas phase: the bare pair is neutral', fontsize=7.5, fontweight='bold', loc='left')
labels = [('B3LYP', 'b3lyp'), ('ωB97X-D3(BJ)', 'wb97xd3bj'), ('M06-2X', 'm062x'), ('PBE0', 'pbe0'),
          ('CAM-B3LYP', 'camb3lyp'), ('DF-MP2', 'dfmp2'), ('SCS-MP2', 'scsmp2')]
ys = list(range(len(labels) + 1))[::-1]
for (lab, key), y in zip(labels, ys):
    v = pcm[key]; bx.barh(y, v, color='#2563eb' if v > 0 else '#7c3aed', height=0.6)
bx.barh(ys[-1], cc, color='#111827', height=0.6)
bx.errorbar(cc, ys[-1], xerr=[[cc - lo], [hi - cc]], color='#9ca3af', capsize=2, lw=0.8)
bx.set_yticks(ys); bx.set_yticklabels([l for l, _ in labels] + ['MP2 + CCSD(T)'], fontsize=6.3)
bx.axvline(0, color='k', lw=0.6); bx.set_xlim(-3.5, 5)
bx.set_xlabel('E(neutral) − E(salt) (kcal mol$^{-1}$)')
bx.set_title('Water (IEF-PCM): the method decides', fontsize=7.5, fontweight='bold', loc='right')
bx.text(4.9, ys[-1], 'neutral lower', ha='right', va='center', fontsize=6)
bx.text(4.9, ys[1] - 0.5, 'salt lower', ha='right', va='center', fontsize=6, color='#2563eb')
fig.suptitle('Benzotriazolium picrate ion pair at low ΔpKa', x=0.5, y=0.985, fontsize=8.5, fontweight='bold')
for ext in ('png', 'pdf', 'tif'):
    fig.savefig(OUT / f'graphical_abstract.{ext}', dpi=600)
print('wrote', OUT)
