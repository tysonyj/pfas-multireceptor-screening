"""
00_verify_ligand_structures.py

Check every ligand SMILES against the registered molecular formula and
molecular weight. Without this check, a structure missing fluorines still
yields physically plausible docking scores while describing a different
molecule entirely.

Usage:  python3 00_verify_ligand_structures.py
"""
import json, os, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D, Descriptors

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lig  = json.load(open(f'{BASE}/config/ligands_anionic.json'))

print(f"{'Compound':<12}{'Formula':>16}{'MW':>9}{'nF':>5}{'charge':>8}{'e-':>7}{'parity':>9}")
print("-" * 68)

bad = []
for name in sorted(lig):
    d = lig[name]
    m = Chem.MolFromSmiles(d['smiles'])
    if m is None:
        bad.append((name, 'SMILES parse failed')); continue
    mh = Chem.AddHs(m)
    formula = D.CalcMolFormula(m)
    mw      = Descriptors.MolWt(m)
    nF      = sum(1 for a in m.GetAtoms() if a.GetSymbol() == 'F')
    q       = d['charge']
    nelec   = sum(a.GetAtomicNum() for a in mh.GetAtoms()) - q

    ok_f  = formula == d['formula']
    ok_mw = abs(mw - d['mw']) < 0.6
    ok_e  = nelec % 2 == 0          # odd electron count = charge/structure mismatch
    if not (ok_f and ok_mw and ok_e):
        bad.append((name, f"formula={ok_f} mw={ok_mw} even_electrons={ok_e}"))

    print(f"{name:<12}{formula:>16}{mw:>9.2f}{nF:>5}{q:>8}{nelec:>7}"
          f"{'even' if ok_e else 'ODD':>9}")

print("-" * 68)
if bad:
    print(f"FAILED: {len(bad)}")
    for n, r in bad: print(f"  {n}: {r}")
    sys.exit(1)
print(f"All {len(lig)} ligands verified (formula, molecular weight, electron parity).")
