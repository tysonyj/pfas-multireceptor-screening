#!/bin/bash
# build_system.sh — build a single GROMACS system
# Usage: bash build_system.sh PXR_PFECHS
set -e

SYS="$1"
[ -z "$SYS" ] && { echo "Usage: bash build_system.sh <system_name>"; exit 1; }

BASE="$(cd "$(dirname "$0")/.." && pwd)"
D="$BASE/systems/$SYS"
MDP="$BASE/mdp"
[ -d "$D" ] || { echo "System directory not found: $D"; exit 1; }

cd "$D"
echo "=============================================="
echo "  Building system $SYS"
echo "=============================================="

# ---------- 1. Ligand parameterization (GAFF2 / acpype) ----------
if [ ! -f LIG.acpype/LIG_GMX.itp ]; then
    echo "[1/6] Ligand parameterization (acpype, GAFF2)"
    # Net charge -1 (sulfonate/carboxylate deprotonated at pH 7.4)
    acpype -i ligand.mol2 -b LIG -n -1 -a gaff2 -c bcc
else
    echo "[1/6] Ligand parameters already present, skipping"
fi

# ---------- 2. Protein topology ----------
echo "[2/6] Protein topology (AMBER99SB-ILDN / TIP3P)"
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -water tip3p -ff amber99sb-ildn -ignh

# ---------- 3. Merge the complex ----------
echo "[3/6] Merging protein and ligand"
python3 "$BASE/scripts/merge_complex.py" \
    --protein protein_processed.gro \
    --ligand  LIG.acpype/LIG_GMX.gro \
    --out     complex.gro
python3 "$BASE/scripts/patch_topol.py" \
    --topol topol.top \
    --lig-itp LIG.acpype/LIG_GMX.itp \
    --lig-name LIG

# ---------- 4. Box + solvation ----------
echo "[4/6] Box construction and solvation"
gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt dodecahedron
gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# ---------- 5. Ions (0.15 M NaCl + neutralization) ----------
echo "[5/6] Adding ions"
gmx grompp -f "$MDP/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 3
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15

# ---------- 6. Energy minimization ----------
echo "[6/6] Energy minimization"
gmx grompp -f "$MDP/em.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 3
gmx mdrun -deffnm em -v

# Index groups for temperature coupling
printf "1 | 13\nname 20 Protein_LIG\n14 | 15 | 16\nname 21 Water_and_ions\nq\n" \
    | gmx make_ndx -f em.gro -o index.ndx || \
    echo "Note: check that the make_ndx group numbers match this system"

echo ""
echo "Done. Next: bash run_md.sh $SYS"
