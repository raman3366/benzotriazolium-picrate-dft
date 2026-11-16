"""
Attach the chemistry to the auto-generated case files: proper compound and fragment names, the ELF/LOL
planes worth mapping, per-element basis overrides, and provenance (CCDC number, source paper).
The geometry itself comes from dft/build_cluster.py and is not touched here.
Usage: python scripts/dft/curate_cases.py
"""
import json, pathlib

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # scripts/ (paths.py)
import paths

META = {   # this dataset's compound only
 'btpc': dict(
   title='Benzotriazolium picrate (BTZP)',
   chemical_name='1H-benzotriazol-1-ium 2,4,6-trinitrophenolate',
   formula='C12H8N6O7', ccdc=835950, doi='10.1016/j.ijleo.2012.06.018',
   fragments={'C6H2N3O7': 'picrate anion', 'C6H6N3': 'benzotriazolium cation'},
   # benzotriazolium is a planar bicycle: the triazole and benzo rings share one plane, so only one map is needed
   planes={'picrate_ring':        ['C6', 'C4', 'C2'],
           'benzotriazolium_ring': ['N4#4_666', 'N5#4_666', 'C7#4_666']},
   notes=('Ion pair assembled from the CIF by symmetry: picrate at 1_555, benzotriazolium at 4_666. The two '
          'moieties of the asymmetric unit are NOT hydrogen bonded to each other (their closest approach is '
          'O4...N5 2.99 A, and N5 carries no hydrogen), so the pair had to be rebuilt from a symmetry image. '
          'The resulting N4-H4N...O1 bond (1.775 A, 154.5 deg) is to the picrate phenolate oxygen, the '
          'expected charge-assisted contact in a picrate salt. The cation\'s second N-H (N6-H6N) donates to a '
          'further symmetry-related picrate and is therefore not represented here.')),
}


def curate():
    for mid, m in META.items():
        path = paths.STRUCTURE / f'{mid}_gas.json'      # draft written by build_cluster.py
        c = json.load(open(path))
        L = c['cif_labels']
        idx = {l: i for i, l in enumerate(L)}
        c['title'] = f"{m['title']}, CCDC {m['ccdc']} -- gas phase"
        c['chemical_name'] = m['chemical_name']
        c['formula'] = m['formula']
        c['ccdc'] = m['ccdc']
        c['source_doi'] = m['doi']
        c['fragments'] = {m['fragments'].get(k, k): v for k, v in c['fragments'].items()}
        c['planes'] = {name: [idx[a] for a in atoms] for name, atoms in m['planes'].items()}
        if m.get('basis_overrides'):
            c['basis_overrides'] = m['basis_overrides']
        if m.get('proton'):
            D, H, A = m['proton']
            c['proton'] = {'H': idx[H], 'donor_N': idx[D], 'acceptor_O': idx[A]}
        c['notes'] = m['notes']
        json.dump(c, open(path, 'w'), indent=2)
        pr = c.get('proton', {})
        print(f"[curate] {mid}: {len(L)} atoms, charge {c['charge']}, "
              f"{len(c['fragments'])} fragments, {len(c['planes'])} ELF/LOL planes, "
              f"basis overrides {c.get('basis_overrides', {})}")
        print(f"          fragments: {', '.join(c['fragments'])}")
        print(f"          planes:    {', '.join(c['planes'])}")
        print(f"          hbonds:    {', '.join(c['hbonds'])}")
        if pr:
            print(f"          proton watched: {L[pr['donor_N']]}-{L[pr['H']]}...{L[pr['acceptor_O']]}")


if __name__ == '__main__':
    curate()
