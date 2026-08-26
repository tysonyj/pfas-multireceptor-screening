#!/bin/bash
# analyze.sh (v2) — trajectory analysis for all four systems plus
#                   protein-ligand interaction energies
#
# Changes from v1:
#   - handles the part files produced by mdrun -noappend (md_prod.part0001.xtc etc.)
#   - generates the final .gro from the last trajectory frame when it is missing
#   - detects the ligand residue name automatically
#   - continues to the next system if any step fails

set -u
BASE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$BASE/analysis"; mkdir -p "$OUT"

for S in PXR_PFECHS PPARa_PFECHS PXR_F-53B PPARa_F-53B; do
  D="$BASE/systems/$S"
  echo ""
  echo "============================================================"
  echo "  $S"
  echo "============================================================"
  [ -d "$D" ] || { echo "  directory not found, skipping"; continue; }
  cd "$D"
  mkdir -p "$OUT/$S"

  # ---------- 1. Locate the trajectory (concatenating part files if needed) ----------
  PARTS=$(ls md_prod.part*.xtc 2>/dev/null | sort)
  NPART=$(echo "$PARTS" | grep -c . || true)
  if [ "$NPART" -gt 1 ]; then
      echo "  concatenating $NPART part files"
      [ -f traj_all.xtc ] || gmx trjcat -f $PARTS -o traj_all.xtc -settime <<< "$(yes c | head -n $NPART)" 2>/dev/null \
                          || gmx trjcat -f $PARTS -o traj_all.xtc -cat
      XTC=traj_all.xtc
  elif [ "$NPART" -eq 1 ]; then
      XTC=$(echo "$PARTS")
  else
      XTC=$(ls md_prod.xtc 2>/dev/null | head -1)
  fi
  [ -n "${XTC:-}" ] && [ -f "$XTC" ] || { echo "  no trajectory found, skipping"; continue; }
  echo "  trajectory: $XTC"

  # ---------- 2. Locate the tpr ----------
  TPR=$(ls md_prod.tpr md_prod.part*.tpr 2>/dev/null | head -1)
  [ -n "${TPR:-}" ] || { echo "  no tpr found, skipping"; continue; }

  # ---------- 3. Locate the final structure (.gro) ----------
  GRO=$(ls md_prod.gro md_prod.part*.gro 2>/dev/null | head -1)
  if [ -z "${GRO:-}" ]; then
      echo "  no final .gro -> extracting the last trajectory frame"
      printf "0\n" | gmx trjconv -s "$TPR" -f "$XTC" -o md_final.gro -dump 999999999 >/dev/null 2>&1
      GRO=md_final.gro
  fi
  echo "  structure: $GRO"

  # ---------- 4. PBC correction ----------
  [ -f nojump.xtc ] || printf "1\n0\n" | \
      gmx trjconv -s "$TPR" -f "$XTC" -o nojump.xtc -pbc mol -center >/dev/null 2>&1

  # ---------- 5. Ligand group number / residue name ----------
  LIGRES=$(awk '$1=="ATOM"||$1=="HETATM"{print substr($0,18,3)}' ligand.pdb 2>/dev/null | sort -u | head -1)
  LIGRES=${LIGRES:-UNL}
  LIGGRP=$(printf "q\n" | gmx make_ndx -f "$GRO" -o /tmp/_probe.ndx 2>&1 | \
           grep -iE "^ *[0-9]+ +$LIGRES" | head -1 | awk '{print $1}')
  LIGGRP=${LIGGRP:-13}
  echo "  ligand: residue $LIGRES, group $LIGGRP"

  # ---------- 6. Metrics ----------
  printf "%s\n%s\n" "$LIGGRP" "$LIGGRP" | gmx rms -s "$TPR" -f nojump.xtc \
      -o "$OUT/$S/rmsd_ligand.xvg" -tu ns >/dev/null 2>&1
  printf "4\n4\n" | gmx rms -s "$TPR" -f nojump.xtc \
      -o "$OUT/$S/rmsd_protein.xvg" -tu ns >/dev/null 2>&1
  printf "1\n%s\n" "$LIGGRP" | gmx mindist -s "$TPR" -f nojump.xtc \
      -od "$OUT/$S/mindist.xvg" -on "$OUT/$S/contacts.xvg" -d 0.4 -tu ns >/dev/null 2>&1
  printf "1\n" | gmx gyrate -s "$TPR" -f nojump.xtc \
      -o "$OUT/$S/gyrate.xvg" >/dev/null 2>&1

  # ---------- 7. Protein-ligand interaction energy ----------
  echo "  interaction energy (70-100 ns)"
  sed "s/^energygrps.*/energygrps      = Protein $LIGRES/" "$BASE/mdp/lie_rerun.mdp" > ie.mdp
  if gmx grompp -f ie.mdp -c "$GRO" -p topol.top -o ie.tpr -maxwarn 5 >/dev/null 2>&1; then
      gmx mdrun -s ie.tpr -rerun nojump.xtc -deffnm ie -nb cpu -ntmpi 1 >/dev/null 2>&1
      if printf "Coul-SR:Protein-%s\nLJ-SR:Protein-%s\n\n" "$LIGRES" "$LIGRES" | \
         gmx energy -f ie.edr -o "$OUT/$S/interaction_energy.xvg" -b 70000 -e 100000 >/dev/null 2>&1; then
          echo "  OK  interaction energy extracted"
      else
          echo "  WARNING  could not extract the energy terms. Available terms:"
          echo | gmx energy -f ie.edr 2>&1 | grep -iE "coul-sr|lj-sr" | head -5
      fi
  else
      echo "  WARNING  grompp failed (could not build ie.tpr)"
  fi

  cd "$BASE"
done

echo ""
echo "============================================================"
echo "  Analysis complete -> $OUT"
echo "  Next: python3 $BASE/scripts/summarize.py"
echo "============================================================"
