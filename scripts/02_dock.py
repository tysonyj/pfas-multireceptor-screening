"""
Production docking (Vina CLI version) — uses the vina executable only, no
Python bindings.

Usage:
    python3 02_dock.py              # when vina is on PATH
    VINA_BIN=/path/to/vina python3 02_dock.py

The run is resumable: restarting after an interruption picks up where it left
off, based on results/docking_results.json.
"""
import json, os, re, csv, time, shutil, subprocess, sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG  = json.load(open(f'{BASE}/config/config_verified.json'))
OUT  = f'{BASE}/results'
os.makedirs(f'{OUT}/poses', exist_ok=True)

EXHAUSTIVENESS = 16
NUM_MODES      = 9
SEED           = 42
CPU            = os.cpu_count() or 4

VINA = os.environ.get('VINA_BIN') or shutil.which('vina')
if not VINA:
    sys.exit(
        "vina executable not found.\n"
        "  1) check with:  which vina\n"
        "  2) if missing, download it:\n"
        "     wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64\n"
        "     chmod +x vina_1.2.5_linux_x86_64\n"
        "     VINA_BIN=./vina_1.2.5_linux_x86_64 python3 02_dock.py"
    )
print(f'vina  : {VINA}')
print(f'cpu   : {CPU}')

# ---- use only targets that passed the validation gate ----
vf = f'{OUT}/redock_validation.json'
if os.path.exists(vf):
    val = json.load(open(vf))
    targets = [t for t in CFG if val.get(t, {}).get('pass')]
    skipped = [t for t in CFG if t not in targets]
    if skipped:
        print(f'excluded (failed validation): {", ".join(skipped)}')
else:
    print('Warning: no validation results file found. Using all targets.')
    targets = list(CFG)

ligands = sorted(f[:-6] for f in os.listdir(f'{BASE}/ligands') if f.endswith('.pdbqt'))
total   = len(ligands) * len(targets)
print(f'{len(targets)} targets x {len(ligands)} ligands = {total} pairs\n')

rf  = f'{OUT}/docking_results.json'
res = json.load(open(rf)) if os.path.exists(rf) else {}
print(f'resuming from {len(res)} existing pairs\n')

AFF = re.compile(r'^\s*1\s+(-?\d+\.\d+)')

def best_affinity(stdout, out_pdbqt):
    for line in stdout.splitlines():
        m = AFF.match(line)
        if m:
            return float(m.group(1))
    if os.path.exists(out_pdbqt):                       # fallback: parse the output file
        for line in open(out_pdbqt):
            if line.startswith('REMARK VINA RESULT'):
                return float(line.split()[3])
    return None

t_start = time.time()
for t in targets:
    c = CFG[t]
    rec = f'{BASE}/receptors/{t}_receptor.pdbqt'
    if not os.path.exists(rec):
        print(f'{t}: receptor file not found, skipping'); continue
    for lg in ligands:
        key = f'{t}|{lg}'
        if key in res:
            continue
        out = f'{OUT}/poses/{t}_{lg}.pdbqt'
        cmd = [VINA,
               '--receptor', rec,
               '--ligand',   f'{BASE}/ligands/{lg}.pdbqt',
               '--center_x', str(c['center'][0]),
               '--center_y', str(c['center'][1]),
               '--center_z', str(c['center'][2]),
               '--size_x',   str(c['size'][0]),
               '--size_y',   str(c['size'][1]),
               '--size_z',   str(c['size'][2]),
               '--exhaustiveness', str(EXHAUSTIVENESS),
               '--num_modes',      str(NUM_MODES),
               '--seed',           str(SEED),
               '--cpu',            str(CPU),
               '--out', out]
        t0 = time.time()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            dG = best_affinity(p.stdout, out)
            if dG is None:
                print(f'{t:<8} {lg:<10} PARSE FAILED: {p.stderr.strip()[:70]}', flush=True)
                res[key] = None
            else:
                res[key] = round(dG, 3)
                n = len(res)
                el = time.time() - t_start
                eta = el / max(n - (total - len(ligands) * len(targets)), 1)
                print(f'{t:<8} {lg:<10} dG={dG:7.3f}  ({time.time()-t0:.0f}s)  [{n}/{total}]', flush=True)
        except subprocess.TimeoutExpired:
            print(f'{t:<8} {lg:<10} TIMEOUT', flush=True); res[key] = None
        except Exception as e:
            print(f'{t:<8} {lg:<10} ERROR {str(e)[:60]}', flush=True); res[key] = None
        json.dump(res, open(rf, 'w'), indent=2)

# ---- CSV matrix ----
with open(f'{OUT}/dG_matrix.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['pfas'] + targets)
    for lg in ligands:
        w.writerow([lg] + [res.get(f'{t}|{lg}', '') for t in targets])

ok = sum(1 for v in res.values() if v is not None)
print(f'\nDone: {ok}/{total} succeeded')
print(f'→ {OUT}/dG_matrix.csv')
