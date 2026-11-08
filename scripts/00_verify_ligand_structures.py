"""
00_verify_ligand_structures.py

Check every ligand SMILES against the CAS-registered molecular formula and
molecular weight. Without this check, a structure missing fluorines still
yields physically plausible docking scores while describing a different
molecule entirely.

The reference values (`cas`, `cas_formula`, `cas_mw` in
config/ligands_anionic.json) are transcribed from the CAS registry entry for
each compound and are deliberately held separately from the SMILES, so that a
SMILES and a reference transcribed from the same erroneous source cannot agree
by construction. A mismatch between the RDKit-derived formula/MW and the
registry values fails the compound.

Usage:  python3 00_verify_ligand_structures.py
"""
import json, os, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as D, Descriptors

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lig  = json.load(open(f'{BASE}/config/ligands_anionic.json'))

print(f"{'Compound':<10}{'CAS':>14}{'Formula (RDKit)':>17}{'CAS formula':>15}"
      f"{'MW':>8}{'CAS MW':>9}{'nF':>4}{'e-':>6}{'result':>9}")
print("-" * 92)

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

    # compare against the CAS registry values, not against the fields that
    # travel with the SMILES
    ok_f  = formula == d['cas_formula']
    ok_mw = abs(mw - d['cas_mw']) < 0.6
    ok_e  = nelec % 2 == 0          # odd electron count = charge/structure mismatch
    ok_sc = formula == d['formula'] # internal self-consistency, reported separately
    if not (ok_f and ok_mw and ok_e):
        bad.append((name, f"cas_formula={ok_f} cas_mw={ok_mw} even_electrons={ok_e}"))

    print(f"{name:<10}{d['cas']:>14}{formula:>17}{d['cas_formula']:>15}"
          f"{mw:>8.2f}{d['cas_mw']:>9.2f}{nF:>4}{nelec:>6}"
          f"{('PASS' if (ok_f and ok_mw and ok_e) else 'FAIL'):>9}")

print("-" * 92)
if bad:
    print(f"FAILED: {len(bad)}")
    for n, r in bad: print(f"  {n}: {r}")
    sys.exit(1)
print(f"All {len(lig)} ligands verified against CAS registry formula and molecular weight,\n"
      f"with even electron count confirmed for every assigned charge state.")
