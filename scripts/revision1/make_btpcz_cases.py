"""
Case files for the second BTPC polymorph (revision item C0): Zeng, Li & Wang, Acta Cryst. E67, o1464 (2011),
doi:10.1107/S1600536811018307, COD 2230422, P2_1/n (the published paper models CCDC 835950, P2_1/c).

Step 1 (geometry, unchanged project tool; the CIF loop lists the defining hydrogen bond with symmetry code '.'):
  python -I dft/build_cluster.py molecules/btpcz/input/2230422.cif btpcz --require-hbond N6:O1 \
      --title "Benzotriazolium picrate, P2_1/n polymorph (Zeng et al. 2011, COD 2230422)"
Step 2 (this script): attach the chemistry exactly as dft/curate_cases.py does for btpc (same fragment names, the
same two ring planes, no basis overrides: 6-311++G(d,p) exists for H, C, N, O) and write the IEF-PCM (water) twin
exactly as dft/pipeline.py make_pcm_case would.  Idempotent: rerunning on its own output gives the same files.
  python dft/revision/make_btpcz_cases.py

Why N6-H6...O1 and not the bifurcated N4-H4...O1/O2: it is the shortest and most linear N-H...O(phenolate) contact
(H...A 1.87 A, N...O 2.770 A, 160 deg), the one the CIF lists within the asymmetric unit, and the analogue of the
single charge-assisted N4-H4N...O1 bond that defines the btpc pair (1.775 A, 2.623 A, 154.5 deg).  Unconstrained,
build_cluster.py prefers the 2_646 image whose N4-H4 bridges O1 and the ortho-nitro O2 (H...A 2.035/2.161 A,
N...O 2.799/2.791 A) because it counts two contacts; that alternative is recorded in the case notes.
"""
import json, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
MID = 'btpcz'
CASES = ROOT / 'molecules' / MID / 'cases'

META = dict(
    title='Benzotriazolium picrate, P2_1/n polymorph (Zeng et al. 2011, COD 2230422)',
    chemical_name='1H-benzotriazol-1-ium 2,4,6-trinitrophenolate',
    formula='C12H8N6O7', cod=2230422, doi='10.1107/S1600536811018307',
    fragments={'C6H2N3O7': 'picrate anion', 'C6H6N3': 'benzotriazolium cation'},
    # same ring planes as btpc (curate_cases.py): picrate C6/C4/C2, benzotriazolium N4/N5/C7
    planes={'picrate_ring': ['C6', 'C4', 'C2'], 'benzotriazolium_ring': ['N4', 'N5', 'C7']},
    proton=('N6', 'H6', 'O1'),
    notes=('Second polymorph of BTPC (revision item C0): Zeng, Li & Wang, Acta Cryst. E67, o1464 (2011), '
           'doi:10.1107/S1600536811018307, COD 2230422, space group P2_1/n, a 14.4113, b 3.7608, c 24.941 A, '
           'beta 90.598 deg, Z 4, 298 K, R1 0.0616. The published paper models the P2_1/c form (CCDC 835950, '
           'case btpc). Here the asymmetric unit IS the hydrogen-bonded ion pair: N6-H6...O1 (phenolate), '
           'H...A 1.870 A, N...O 2.770 A, 160.0 deg, matching the deposited _geom_hbond loop (2.770(3) A, '
           '160(3) deg). The cation carries H on N4 and N6 and none on N5, as in btpc (1H,3H-benzotriazolium); '
           'the deposited N-H distances are refined X-ray values (N4-H4 0.85, N6-H6 0.94 A). The cation\'s '
           'other N-H (N4-H4) donates a bifurcated bond to O1 and the ortho-nitro O2 of a 2_656 picrate '
           '(N...O 2.799 and 2.791 A); it is not represented here. Unconstrained, build_cluster.py would have '
           'picked that bifurcated image (two contacts outscore one), so the pair was built with '
           '--require-hbond N6:O1. No CCDC number is given in the COD CIF.'),
)


def curate(c):
    L = c['cif_labels']
    idx = {l: i for i, l in enumerate(L)}
    c['title'] = f"{META['title']} -- gas phase"
    c['chemical_name'] = META['chemical_name']
    c['formula'] = META['formula']
    c.pop('ccdc', None)
    c['cod'] = META['cod']
    c['source_doi'] = META['doi']
    c['fragments'] = {META['fragments'].get(k, k): v for k, v in c['fragments'].items()}
    c['planes'] = {name: [idx[a] for a in atoms] for name, atoms in META['planes'].items()}
    D, H, A = META['proton']
    c['proton'] = {'H': idx[H], 'donor_N': idx[D], 'acceptor_O': idx[A]}
    c['notes'] = META['notes']
    c.pop('solvent', None)
    return c


def pcm_twin(gas):
    """dft/pipeline.py make_pcm_case: same content, solvent water, title and notes adjusted."""
    new = json.loads(json.dumps(gas))
    new['solvent'] = 'water'
    t = new['title']
    new['title'] = t.replace('gas phase', 'IEF-PCM (water)') if 'gas phase' in t else t + ' -- IEF-PCM (water)'
    new['notes'] = (new.get('notes', '') + ' | Rerun with an IEF-PCM (water) continuum because the continuum '
                    'model is the primary result for an ionic salt.').strip()
    return new


def check(c, nat):
    assert len(c['cif_labels']) == nat, 'labels vs atoms'
    allf = sorted(i for v in c['fragments'].values() for i in v)
    assert allf == list(range(nat)), 'fragments must partition the atoms'
    assert set(c['fragments']) == {'picrate anion', 'benzotriazolium cation'}, c['fragments'].keys()
    assert c['charge'] == 0 and c['spin'] == 0
    p = c['proton']
    L = c['cif_labels']
    assert (L[p['donor_N']], L[p['H']], L[p['acceptor_O']]) == META['proton']
    for name, (d, h, a) in c['hbonds'].items():
        assert max(d, h, a) < nat
    return True


if __name__ == '__main__':
    src = CASES / f'{MID}_gas.json'
    gas = curate(json.load(open(src)))
    xyz = ROOT / gas['xyz']
    nat = int(open(xyz).readline())
    check(gas, nat)
    pcm = pcm_twin(gas)
    check(pcm, nat)
    json.dump(gas, open(src, 'w'), indent=2)
    json.dump(pcm, open(CASES / f'{MID}_pcm.json', 'w'), indent=2)
    L = gas['cif_labels']; p = gas['proton']
    print(f"[btpcz] {nat} atoms; fragments {', '.join(f'{k} ({len(v)})' for k, v in gas['fragments'].items())}")
    print(f"[btpcz] proton watched: {L[p['donor_N']]}-{L[p['H']]}...{L[p['acceptor_O']]}; hbonds {list(gas['hbonds'])}")
    print(f"[btpcz] planes {gas['planes']}; wrote {src.relative_to(ROOT)} and {MID}_pcm.json")
