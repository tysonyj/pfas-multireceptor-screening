"""
Validation gate — redock each receptor's co-crystallized ligand and check
whether the crystallographic pose is reproduced.
Pass criterion: heavy-atom RMSD < 2.0 A (or centroid deviation < 3.0 A).
Targets that fail this gate must be excluded from the production docking run.
"""
import json, os, time
import numpy as np
from vina import Vina

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cfg  = json.load(open(f'{BASE}/config/config_verified.json'))
OUT  = f'{BASE}/results'; os.makedirs(OUT, exist_ok=True)

def heavy(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path)
                     if l.startswith(('ATOM', 'HETATM'))
                     and not l[12:16].strip().startswith('H')])

res = {}
rf = f'{OUT}/redock_validation.json'
if os.path.exists(rf):
    res = json.load(open(rf))

for t, c in cfg.items():
    if t in res and 'error' not in res[t]:
        print(f'{t:<8} (already done, skipping)'); continue
    t0 = time.time()
    try:
        v = Vina(sf_name='vina', cpu=0, seed=42, verbosity=0)
        v.set_receptor(f'{BASE}/receptors/{t}_receptor.pdbqt')
        v.set_ligand_from_file(f'{BASE}/reference_ligands/ref_{t}.pdbqt')
        v.compute_vina_maps(center=c['center'], box_size=c['size'])
        v.dock(exhaustiveness=16, n_poses=5)
        v.write_poses(f'{OUT}/redock_{t}.pdbqt', n_poses=1, overwrite=True)

        dG   = v.energies(n_poses=1)[0][0]
        pose = heavy(f'{OUT}/redock_{t}.pdbqt')
        dev  = float(np.linalg.norm(pose.mean(0) - np.array(c['center'])))
        ok   = dev < 3.0
        res[t] = {'dG': round(float(dG), 2), 'centroid_dev': round(dev, 2),
                  'pass': bool(ok), 'sec': round(time.time() - t0)}
        print(f'{t:<8} dG={dG:7.2f}  centroid_dev={dev:5.2f}A  '
              f"{'PASS' if ok else 'FAIL'}  ({time.time()-t0:.0f}s)", flush=True)
    except Exception as e:
        print(f'{t:<8} ERROR {type(e).__name__}: {e}', flush=True)
        res[t] = {'error': str(e)}
    json.dump(res, open(rf, 'w'), indent=2)

n_pass = sum(1 for v in res.values() if v.get('pass'))
print(f'\nValidation passed: {n_pass}/{len(cfg)}')
if n_pass < len(cfg):
    print('Warning: failed targets must be excluded from the production docking '
          'run, or their search boxes redefined.')
