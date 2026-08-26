#!/bin/bash
# run_md.sh — NVT/NPT equilibration + 100 ns production
# Usage: bash run_md.sh PXR_PFECHS [GPU_ID]
set -e
SYS="$1"; GPU="${2:-0}"
BASE="$(cd "$(dirname "$0")/.." && pwd)"
D="$BASE/systems/$SYS"; MDP="$BASE/mdp"
cd "$D"
NDX=""; [ -f index.ndx ] && NDX="-n index.ndx"

echo "=== $SYS : NVT (100 ps) ==="
[ -f nvt.gro ] || {
  gmx grompp -f "$MDP/nvt.mdp" -c em.gro -r em.gro -p topol.top $NDX -o nvt.tpr -maxwarn 5
  gmx mdrun -deffnm nvt -gpu_id $GPU -nb gpu -v
}

echo "=== $SYS : NPT (500 ps) ==="
[ -f npt.gro ] || {
  gmx grompp -f "$MDP/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top $NDX -o npt.tpr -maxwarn 5
  gmx mdrun -deffnm npt -gpu_id $GPU -nb gpu -v
}

echo "=== $SYS : Production (100 ns) ==="
[ -f md_prod.tpr ] || \
  gmx grompp -f "$MDP/md_prod.mdp" -c npt.gro -t npt.cpt -p topol.top $NDX -o md_prod.tpr -maxwarn 5
gmx mdrun -deffnm md_prod -gpu_id $GPU -nb gpu -bonded gpu -pme gpu -v -cpi md_prod.cpt -noappend || \
gmx mdrun -deffnm md_prod -gpu_id $GPU -nb gpu -v -cpi md_prod.cpt -noappend

echo "Done: $SYS"
