#!/bin/bash
# rebuild_pxr.sh (v2) — rebuild the PXR systems
#
# Problem with v1: protein.pdb is renumbered when Open Babel converts PDBQT to
#                  PDB, so cutting at the crystallographic residue number (198)
#                  removes the wrong part of the chain.
#                  -> detect the fragments directly from backbone connectivity
#                     instead of relying on residue numbering.
#
# Usage:  bash rebuild_pxr.sh          (prepare + verify)
#         bash rebuild_pxr.sh --run    (prepare + run MD)

set -e
BASE="$(cd "$(dirname "$0")/.." && pwd)"
RUN=${1:-}

for SYS in PXR_PFECHS PXR_F-53B; do
  OLD="$BASE/systems/$SYS"
  NEW="$BASE/systems/${SYS}_v2"
  echo ""
  echo "============================================================"
  echo "  $SYS  →  ${SYS}_v2"
  echo "============================================================"
  [ -d "$OLD" ] || { echo "  source directory not found"; continue; }
  mkdir -p "$NEW"

  python3 - "$OLD" "$NEW" << 'PY'
import sys, numpy as np
from collections import OrderedDict
old, new = sys.argv[1], sys.argv[2]

# ---------- Collect atoms by residue ----------
res = OrderedDict()
for l in open(f'{old}/protein.pdb'):
    if l.startswith('ATOM'):
        key = (l[21], l[22:27])                # chain + resSeq+icode
        res.setdefault(key, []).append(l)
keys = list(res)
print(f"  input: {len(keys)} residues / {sum(len(v) for v in res.values())} atoms")

# ---------- Split into fragments using consecutive C(i)-N(i+1) distances ----------
def atom(lines, name):
    for l in lines:
        if l[12:16].strip() == name:
            return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return None

frag, cur = [], [keys[0]]
for a, b in zip(keys, keys[1:]):
    C, N = atom(res[a], 'C'), atom(res[b], 'N')
    d = np.linalg.norm(C - N) if (C is not None and N is not None) else 99.0
    if d > 2.0:                                 # peptide bond is ~1.33 A
        frag.append(cur); cur = [b]
    else:
        cur.append(b)
frag.append(cur)
frag.sort(key=len, reverse=True)
print(f"  {len(frag)} fragments: " + ", ".join(f"{len(f)} residues" for f in frag))

if len(frag) == 1:
    print("  no chain break — rebuild not required")
    sys.exit(0)

keep = set(frag[0])                             # keep only the largest fragment
with open(f'{new}/protein.pdb', 'w') as fh:
    for k in keys:
        if k in keep:
            fh.writelines(res[k])
    fh.write("END\n")
print(f"  kept: {len(frag[0])} residues / {sum(len(res[k]) for k in keep)} atoms")
print(f"  removed: " + ", ".join(f"{len(f)} residues" for f in frag[1:]))

# ---------- Copy the ligand and verify the binding site ----------
import shutil
for f in ('ligand.pdb', 'ligand.mol2', 'info.json'):
    try: shutil.copy(f'{old}/{f}', f'{new}/{f}')
    except FileNotFoundError: pass

def coords(p, recs):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(p) if l.startswith(recs)])
P = coords(f'{new}/protein.pdb', 'ATOM')
L = coords(f'{new}/ligand.pdb', ('ATOM', 'HETATM'))
D = np.linalg.norm(P[:, None] - L[None], axis=2)
n45 = int((D < 4.5).sum())
print(f"  check: minimum ligand-protein distance {D.min():.2f} A, {n45} contacts within 4.5 A")
print("  " + ("OK  binding site preserved" if (D.min() < 4.0 and n45 > 15)
               else "WARNING  too few contacts — inspect before proceeding"))
PY
done

echo ""
echo "============================================================"
if [ "$RUN" == "--run" ]; then
  for SYS in PXR_PFECHS_v2 PXR_F-53B_v2; do
    echo ""; echo "###  $SYS  ###"
    bash "$BASE/scripts/build_system.sh" "$SYS"
    bash "$BASE/scripts/run_md.sh" "$SYS" 0
  done
else
  echo "  If the check passes, run:  bash rebuild_pxr.sh --run"
fi
echo "============================================================"
