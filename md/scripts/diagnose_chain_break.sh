#!/bin/bash
# diagnose_chain_break.sh — determine whether the anomalous backbone RMSD of the
# PXR systems is a reference-structure artefact or a genuine structural change.
#
# Background. The PXR crystal structure 1ILH is missing residues 178-197, so
# GROMACS treats the protein as two separate molecules. Two observations rule
# out physical dissociation:
#   - the two segments (142-177 / 198-431) stay in contact throughout, with a
#     minimum separation of 0.19 nm
#   - the radius of gyration of PXR_PFECHS drifts by only +2.3%
# Yet the backbone RMSD measured against the tpr is 2.86 nm from t = 0.
#
# The test. Recompute the RMSD against the first trajectory frame instead of the
# tpr reference coordinates:
#   RMSD vs frame 0 < 0.5 nm  -> tpr reference-coordinate artefact; the
#                                simulation itself is fine
#   RMSD vs frame 0 still large -> a genuine large-scale structural change
#
# Usage: bash diagnose_chain_break.sh
# Output: analysis_diag/<system>/

set -u
BASE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$BASE/analysis_diag"; mkdir -p "$OUT"

for S in PXR_PFECHS PXR_F-53B PPARa_PFECHS PPARa_F-53B; do
  D="$BASE/systems/$S"; [ -d "$D" ] || continue
  cd "$D"; mkdir -p "$OUT/$S"
  echo ""
  echo "=========================  $S  ========================="

  XTC=$(ls md_prod.part*.xtc md_prod.xtc 2>/dev/null | head -1)
  TPR=$(ls md_prod.tpr md_prod.part*.tpr 2>/dev/null | head -1)
  [ -n "${XTC:-}" ] && [ -n "${TPR:-}" ] || { echo "  files not found"; continue; }

  # Subsample to 1 ns + whole + nojump
  printf "0\n" | gmx trjconv -s "$TPR" -f "$XTC" -o _s.xtc -dt 1000 -pbc whole >/dev/null 2>&1
  printf "0\n" | gmx trjconv -s "$TPR" -f _s.xtc -o _sn.xtc -pbc nojump >/dev/null 2>&1

  # Extract the first trajectory frame as an alternative reference structure
  printf "0\n" | gmx trjconv -s "$TPR" -f _sn.xtc -o _ref.gro -dump 0 >/dev/null 2>&1

  # (A) RMSD against the tpr
  printf "4\n4\n" | gmx rms -s "$TPR" -f _sn.xtc -o "$OUT/$S/rmsd_vs_tpr.xvg" \
      -tu ns -fit rot+trans >/dev/null 2>&1
  # (B) RMSD against frame 0   <- the decisive test
  printf "4\n4\n" | gmx rms -s _ref.gro -f _sn.xtc -o "$OUT/$S/rmsd_vs_frame0.xvg" \
      -tu ns -fit rot+trans >/dev/null 2>&1
  # (C) Radius of gyration
  printf "1\n" | gmx gyrate -s "$TPR" -f _sn.xtc -o "$OUT/$S/gyrate_fixed.xvg" >/dev/null 2>&1

  python3 - "$OUT/$S" << 'PY'
import sys, os, numpy as np
d=sys.argv[1]
def load(f):
    p=os.path.join(d,f)
    if not os.path.exists(p): return None
    v=[]
    for l in open(p):
        if l.startswith(('#','@')): continue
        try: v.append([float(x) for x in l.split()])
        except: pass
    return np.array(v) if v else None
a=load('rmsd_vs_tpr.xvg'); b=load('rmsd_vs_frame0.xvg'); g=load('gyrate_fixed.xvg')
if a is not None:
    y=a[:,1]; print(f"    RMSD vs tpr     : start {y[0]:6.3f}  end {y[-10:].mean():6.3f}  max {y.max():6.3f} nm")
if b is not None:
    y=b[:,1]; m=y[-10:].mean()
    print(f"    RMSD vs frame 0 : start {y[0]:6.3f}  end {m:6.3f}  max {y.max():6.3f} nm  <-")
    print(f"    -> {'reference-coordinate artefact; protein is fine' if m<0.5 else 'genuine structural change (>0.5 nm)'}")
if g is not None:
    y=g[:,1]; print(f"    Rg              : start {y[0]:.3f} -> end {y[-10:].mean():.3f} nm "
                    f"({100*(y[-10:].mean()/y[0]-1):+.1f}%)")
PY
  rm -f _s.xtc _sn.xtc _ref.gro
  cd "$BASE"
done
echo ""
echo "======================================================="
echo "  The RMSD against frame 0 is the decisive measurement"
echo "======================================================="
