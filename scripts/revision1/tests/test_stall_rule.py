"""
Calibration of opt_d3.py's stall stop (--negligible-stop) and cap rule on the real btpc_pcm optimisation
(data/r2/molecules/btpc/results/btpc_pcm/logs/optimize.log: 71 evaluations, accepted by hand on forces at step 70).

Checks
  1. the default stall rule (window 10, energy spread < 1e-6 Ha, Gaussian gradient thresholds) never fires on that
     run: the forces were inside the thresholds from step 6, but the energy kept falling (0.66 mHa more by step 70,
     with steps of up to 0.3 A), so stopping there would have been premature;
  2. a forces-only style rule (window 5, spread < 1e-5 Ha, i.e. the proton_scan.py maxsteps rule applied at every
     step) would have stopped at about step 12, 0.66 mHa above the accepted structure -- the reason it is used only
     at a cap;
  3. on a synthetic stall appended to the log (trust radius collapsed: |dE| ~ 5e-8, forces 1e-4) the default rule
     fires after 10 evaluations.
Usage: THREADS=1 python dft/revision/tests/test_stall_rule.py   (exit 0 = all checks pass)
"""
import sys, re, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import rev_common as rc          # noqa: F401  (pins threads before numpy)
from opt_d3 import stalled, cap_rule

LOG = rc.ROOT / 'data' / 'r2' / 'molecules' / 'btpc' / 'results' / 'btpc_pcm' / 'logs' / 'optimize.log'
GRMS, GMAX = 3.0e-4, 4.5e-4


def parse(path):
    hist = []
    for line in open(path, errors='replace'):
        line = re.sub(r'\x1b\[[0-9;]*m', '', line)
        m0 = re.match(r'Step\s+0 : Gradient = ([\d.e+-]+)/([\d.e+-]+) \(rms/max\) Energy = (-?[\d.]+)', line)
        m = re.match(r'Step\s+(\d+) : .*Grad = ([\d.e+-]+)/([\d.e+-]+) \(rms/max\) E \(change\) = (-?[\d.]+)', line)
        if m0:
            hist.append({'step': 0, 'grms': float(m0.group(1)), 'gmax': float(m0.group(2)), 'E_Ha': float(m0.group(3))})
        elif m:
            hist.append({'step': int(m.group(1)), 'grms': float(m.group(2)), 'gmax': float(m.group(3)),
                         'E_Ha': float(m.group(4))})
    return hist


def first(hist, rule):
    for k in range(1, len(hist) + 1):
        if rule(hist[:k]):
            return k - 1
    return None


def main():
    h = parse(LOG)
    ok = True
    print(f'{LOG.relative_to(rc.ROOT)}: {len(h)} evaluations, E {h[0]["E_Ha"]:.8f} -> {h[-1]["E_Ha"]:.8f} Ha')
    inside = [x['step'] for x in h if x['grms'] <= GRMS and x['gmax'] <= GMAX]
    print(f'  gradient criteria met at {len(inside)} of {len(h)} evaluations (first at step {inside[0]})')
    s_default = first(h, lambda x: stalled(x, 10, GRMS, GMAX, 1e-6))
    print(f'  default stall rule (W=10, 1e-6 Ha) fires at: {s_default}  (expected: never)')
    ok &= s_default is None
    s_loose = first(h, lambda x: stalled(x, 5, GRMS, GMAX, 1e-5))
    if s_loose is not None:
        print(f'  forces-only rule (W=5, 1e-5 Ha) would fire at step {s_loose}, E {h[s_loose]["E_Ha"]:.8f} Ha, '
              f'{(h[s_loose]["E_Ha"] - h[-1]["E_Ha"]) * 1e3:.2f} mHa above step {h[-1]["step"]}')
    ok &= s_loose is not None and s_loose < 30
    print(f'  cap rule at the last step ({h[-1]["step"]}): {cap_rule(h, GRMS, GMAX)} '
          f'(energy spread over the last 5: {(max(x["E_Ha"] for x in h[-5:]) - min(x["E_Ha"] for x in h[-5:])):.2e} Ha)')
    # synthetic stall after the last step
    e0 = h[-1]['E_Ha']
    syn = h + [{'step': h[-1]['step'] + 1 + i, 'grms': 1.0e-4, 'gmax': 2.0e-4, 'E_Ha': e0 - 5e-8 * (i % 2)}
               for i in range(12)]
    s_syn = first(syn, lambda x: stalled(x, 10, GRMS, GMAX, 1e-6))
    print(f'  synthetic stall appended at step {h[-1]["step"] + 1}: default rule fires at step {s_syn} (expected '
          f'{h[-1]["step"] + 9} or {h[-1]["step"] + 10}: a window of 10 that may include the last real step)')
    ok &= s_syn is not None and h[-1]['step'] + 9 <= s_syn <= h[-1]['step'] + 10
    print('PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
